import glob
import shutil
import regions

import os
import time
import numpy as np
import warnings
import datetime
from astropy import units as u
from astropy.stats import mad_std
from astropy.table import Table
from astropy import log
from spectral_cube import SpectralCube
from spectral_cube.utils import NoBeamError

from casa_formats_io import Table as casaTable

from pathlib import Path

from aces import conf
basepath = conf.basepath

tbldir = Path(f'{basepath}/reduction_ACES/aces/data/tables')

dataroot = f'{basepath}/data/2021.1.00172.L'

global then
then = time.time()

# number of columns in the table.  This has to be checked globally because any
# edits to the underlying table will result in inconsistencies & crashes.
# i.e., any time we add or remove a column from the table, we need to make sure
# this number specifies the number of columns
NCOLS = 44


def dt(message=""):
    global then
    now = time.time()
    print(f"Elapsed: {now - then:0.1g}.  {message}", flush=True)
    then = now


def main(num_workers=None):

    if os.getenv('NO_PROGRESSBAR') is None and not (os.getenv('ENVIRON') == 'BATCH'):
        from dask.diagnostics import ProgressBar
        pbar = ProgressBar()
        pbar.register()

    print(f"TMPDIR = {os.environ.get('TMPDIR')}")

    # Dask writes some log stuff; let's make sure it gets written on local scratch
    # or in a blue drive and not on orange
    if os.getenv('TMPDIR'):
        os.chdir(os.getenv('TMPDIR'))

    threads = int(os.getenv('DASK_THREADS') or os.getenv('SLURM_NTASKS'))
    print(f"Using {threads} threads.")

    spws = {3: list(range(5)), }

    suffix = '.image'

    dt(f"PID = {os.getpid()}")

    if threads:
        # try dask.distrib again
        from dask.distributed import Client, LocalCluster
        import dask

        mem_mb = int(os.getenv('SLURM_MEM_PER_NODE'))
        print(f"Threads was set to {threads}", flush=True)

        try:
            nthreads = int(threads)
            #memlimit = f'{0.8 * int(mem_mb) / int(nthreads)}MB'
            memlimit = f'{0.4 * int(mem_mb)}MB'
            if nthreads > 1:
                num_workers = nthreads
                scheduler = 'threads'
            else:
                scheduler = 'synchronous'
        except (TypeError, ValueError) as ex:
            print(f"Exception raised when creating scheduler: {ex}", flush=True)
            nthreads = 1
            scheduler = 'synchronous'
    else:
        nthreads = 1
        scheduler = 'synchronous'

    target_chunksize = int(1e8)
    print(f"Target chunk size = {target_chunksize} (log10={np.log10(target_chunksize)})", flush=True)

    print(f"Using scheduler {scheduler} with {nthreads} threads", flush=True)
    time.sleep(1)
    print("Slept for 1s", flush=True)

    cwd = os.getcwd()
    os.chdir(basepath)
    print(f"Changed from {cwd} to {basepath}, now running cube stats assembly", flush=True)

    colnames_apriori = ['Field', 'Config', 'spw', 'suffix', 'filename', 'bmaj', 'bmin', 'bpa', 'wcs_restfreq', 'minfreq', 'maxfreq']
    colnames_fromheader = ['imsize', 'cell', 'threshold', 'niter', 'pblimit', 'pbmask', 'restfreq', 'nchan',
                           'width', 'start', 'chanchunks', 'deconvolver', 'weighting', 'robust', 'git_version', 'git_date', 'version']
    colnames_stats = list(('min max std sum mean'.split() +
                      'min_K max_K std_K sum_K mean_K'.split() +
                      'lowmin lowmax lowstd lowmadstd lowsum lowmean'.split()
                      )) + ['timestamp']

    colnames = colnames_apriori + colnames_fromheader + colnames_stats
    # sanity check to make sure I didn't mis-count things above
    assert len(colnames) == NCOLS

    def try_qty(x):
        try:
            return u.Quantity(x)
        except Exception:
            return list(x)

    def save_tbl(rows, colnames):
        columns = list(map(try_qty, zip(*rows)))
        tbl = Table(columns, names=colnames)
        tbl.write(tbldir / 'feathered_cube_stats.ecsv', overwrite=True)
        tbl.write(tbldir / 'feathered_cube_stats.ipac', format='ascii.ipac', overwrite=True)
        tbl.write(tbldir / 'feathered_cube_stats.html', format='ascii.html', overwrite=True)
        tbl.write(tbldir / 'feathered_cube_stats.tex', overwrite=True)
        tbl.write(tbldir / 'feathered_cube_stats.js.html', format='jsviewer', overwrite=True)
        return tbl

    if os.getenv('START_FROM_CACHED') == 'False':
        start_from_cached = False  # TODO: make a parameter
    else:
        start_from_cached = True
        print(f"Starting from cached file {tbldir / 'feathered_cube_stats.ecsv'}")
    tbl = None
    if start_from_cached and os.path.exists(tbldir / 'feathered_cube_stats.ecsv'):
        tbl = Table.read(tbldir / 'feathered_cube_stats.ecsv')
        print(tbl)

        if np.any(np.isnan(tbl['min'])):
            print(f"There are {np.isnan(tbl['min']).sum()} NaNs in the table.  Will recompute those. len(tbl)={len(tbl)}")
            tbl = tbl[np.isfinite(tbl['min'])]
            print(f"Cut-down table length = {len(tbl)}")

        if len(tbl.colnames) != NCOLS:
            warnings.warn("Cached file is BAD!  Moving it.")
            shutil.move(tbldir / 'feathered_cube_stats.ecsv',
                        tbldir / f'feathered_cube_stats_{datetime.datetime.now().isoformat()}.ecsv')
            rows = []
        else:
            rows = [[tbl[cn].quantity[ii]
                     if tbl[cn].unit not in (None, u.dimensionless_unscaled)
                     else tbl[cn][ii] for cn in tbl.colnames]
                    for ii in range(len(tbl))]
    else:
        rows = []

    cache_stats_file = open(tbldir / "feathered_cube_stats.txt", 'w')

    filelists = (
                 glob.glob(f"{basepath}/upload/HCOp_feather_images/MOPRA_12M_7M_feather/*statcont.contsub.fits") +
                 glob.glob(f"{basepath}/upload/Feather_12m_7m_TP/SPW*/cubes/*statcont.contsub.fits") +
                 glob.glob(f"{basepath}/MOPRA_12M_7M_feather/*statcont.contsub.fits") +
                 glob.glob(f"{basepath}/upload/HCOp_feather_images/*statcont.contsub.fits")
                 )
    ufilelist = np.unique(filelists)
    print(f"Found {len(ufilelist)} unique files out of {len(filelists)} total files")

    for fullpath in ufilelist:
        print(f"Processing {fullpath}")
        fn = fullpath
        basename = os.path.basename(fullpath)
        field = basename[9:11]

        config = 'TP_7M_12M_'
        suffix = '.image.statcont.contsub.fits'
        if 'SPW' in fullpath:
            spw = fullpath.split("SPW")[1][:2]
        elif 'hco+' in fullpath:
            spw = '29'
        elif 'hnco' in fullpath:
            spw = '31'
        else:
            print(f"ERROR: Unknown spw for {fullpath}")
            spw = '99'

        if tbl is not None:
            row_matches = ((tbl['Field'] == field) &
                           (tbl['Config'] == config) &
                           (tbl['spw'] == spw) &
                           (tbl['suffix'] == suffix))
            if any(row_matches) and not np.all(~np.isfinite(tbl[row_matches]['min'])):
                print(f"Skipping {fullpath} as complete: {tbl[row_matches]}", flush=True)
                continue

        print(f"Beginning field {field} config {config} spw {spw} suffix {suffix}", flush=True)
        print(f"File: '{fn}'", flush=True)

        if 'fits' in fn:
            cube = SpectralCube.read(fn, format='fits', use_dask=True)
        else:
            cube = SpectralCube.read(fn, format='casa_image', target_chunksize=target_chunksize)

        sched = cube.use_dask_scheduler(scheduler=scheduler, num_workers=num_workers)
        # print(f"Rechunking {cube} to tmp dir", flush=True)
        # cube = cube.rechunk(save_to_tmp_dir=True)
        # cube.use_dask_scheduler(scheduler)

        try:
            if hasattr(cube, 'beam'):
                beam = cube.beam
        except NoBeamError as ex:
            print(f"Beam not found: {ex}")
            continue

        if hasattr(cube, 'beams'):
            beams = cube.beams
            # use the middle-ish beam
            beam = beams[len(beams) // 2]

        print(f"Beam: {beam}, {beam.major}, {beam.minor}", flush=True)

        history = {}
        history['imsize'] = str(cube.shape[1:])
        history['cell'] = str([x.to(u.arcsec).to_string() for x in cube.wcs.celestial.proj_plane_pixel_scales()])
        history['restfreq'] = float(cube.wcs.wcs.restfrq)
        history['nchan'] = int(cube.shape[0])

        with sched:
            # mask to select the channels with little/less emission
            meanspec = cube.mean(axis=(1, 2))
            lowsignal = meanspec < np.nanpercentile(meanspec, 25)

            print(f"Low-signal region selected {lowsignal.sum()} channels out of {lowsignal.size}."
                  f" ({lowsignal.sum() / lowsignal.size * 100:0.2f}) %")

            assert lowsignal.sum() > 0
            assert lowsignal.sum() < lowsignal.size

            noiseest_cube = cube

            dt(cube)
            dt(noiseest_cube)

            minfreq = cube.with_spectral_unit(u.GHz).spectral_axis.min()
            maxfreq = cube.with_spectral_unit(u.GHz).spectral_axis.max()
            restfreq = cube.wcs.wcs.restfrq

            # print("getting filled data")
            # data = cube._get_filled_data(fill=np.nan)
            # print("finished getting filled data")
            # del data

            # try this as an experiment?  Maybe it's statistics that causes problems?
            #print(f"Computing cube mean with scheduler {scheduler} and sched args {cube._scheduler_kwargs}")
            #mean = cube.mean()
            dt(f"Computing cube statistics with scheduler {scheduler} and sched args {cube._scheduler_kwargs}")
            stats = cube.statistics()
            dt("finished cube stats")
            min = stats['min']
            max = stats['max']
            if np.isnan(min) or np.isnan(max):
                raise ValueError("Cube stats reported a NAN min/max")
            std = stats['sigma']
            sum = stats['sum']
            mean = stats['mean']

            faintstats = noiseest_cube.with_mask(lowsignal[:, None, None]).statistics()
            dt("finished low-signal cube stats")
            lowmin = stats['min']
            lowmax = stats['max']
            lowstd = stats['sigma']
            lowsum = stats['sum']
            lowmean = stats['mean']
            dt("Doing low-signal cube mad-std")
            # we got warnings that this was making large chunks.  Not sure there's an alternative here?
            #with dask.config.set(**{'array.slicing.split_large_chunks': True}): # try to split more
            with dask.config.set(**{'array.slicing.split_large_chunks': False}):  # silence warning
                flatdata = noiseest_cube.with_mask(lowsignal[:, None, None]).flattened()
            dt("Loaded flatdata")
            lowmadstd = mad_std(flatdata)
            dt("Done low-signal cube mad-std")

        del stats
        del faintstats

        del cube

        jtok_equiv = beam.jtok_equiv(u.Quantity(minfreq + maxfreq, u.Hz) / 2)

        row = ([field, config, spw, suffix, os.path.basename(fn),
                beam.major.to(u.arcsec), beam.minor.to(u.arcsec), beam.pa,
                u.Quantity(restfreq, u.Hz), u.Quantity(minfreq, u.Hz), u.Quantity(maxfreq, u.Hz)] +
               [history[key] if key in history else '' for key in colnames_fromheader] +
               [min, max, std, sum, mean] +
               list(map(lambda x: u.Quantity(x).to(u.K, jtok_equiv), [min, max, std, sum, mean])) +
               [lowmin, lowmax, lowstd, lowmadstd, lowsum, lowmean, datetime.datetime.now().isoformat()]
               )
        assert len(row) == len(colnames) == NCOLS
        rows.append(row)

        cache_stats_file.write(" ".join(map(str, row)) + "\n")
        cache_stats_file.flush()
        print(f'len(rows): {len(rows)}, len(colnames): {len(colnames)}')
        tbl = save_tbl(rows, colnames)

        if np.any(np.isnan(tbl['min'])):
            print(f"After processing {fn}, there are {np.isnan(tbl['min']).sum()} NaNs in the table.")

    cache_stats_file.close()

    print(tbl)

    os.chdir(cwd)

    if threads and nthreads > 1 and 'client' in locals():
        client.close()  # noqa
