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

from aces.analysis.imstats import get_psf_secondpeak, get_noise_region
from aces.retrieval_scripts.mous_map import get_mous_to_sb_mapping
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
NCOLS = 50


def dt(message=""):
    global then
    now = time.time()
    print(f"Elapsed: {now-then:0.1g}.  {message}", flush=True)
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
            memlimit = f'{0.4*int(mem_mb)}MB'
            if nthreads > 1:
                num_workers = nthreads
                scheduler = 'threads'

            # the cluster approach turned out to be very inefficient
            elif False:
                print(f"nthreads = {nthreads} > 1, so starting a LocalCluster with memory limit {memlimit}", flush=True)
                #scheduler = 'threads'
                # set up cluster and workers
                cluster = LocalCluster(n_workers=1,
                                       threads_per_worker=int(nthreads),
                                       memory_target_fraction=0.60,
                                       memory_spill_fraction=0.65,
                                       memory_pause_fraction=0.7,
                                       # memory_terminate_fraction=0.9,
                                       memory_limit=memlimit,
                                       silence_logs=False,  # https://stackoverflow.com/questions/58014417/seeing-logs-of-dask-workers
                                       )
                print(f"Created a cluster {cluster}", flush=True)
                client = Client(cluster)
                print(f"Created a client {client}", flush=True)
                scheduler = client
                # https://github.com/dask/distributed/issues/3519
                # https://docs.dask.org/en/latest/configuration.html
                dask.config.set({"distributed.workers.memory.terminate": 0.75})
                print(f"Started dask cluster {client} with mem limit {memlimit}", flush=True)
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
    basepath = dataroot
    os.chdir(basepath)
    print(f"Changed from {cwd} to {basepath}, now running cube stats assembly", flush=True)

    colnames_apriori = ['Field', 'Config', 'spw', 'suffix', 'filename', 'bmaj', 'bmin', 'bpa', 'wcs_restfreq', 'minfreq', 'maxfreq']
    colnames_fromheader = ['imsize', 'cell', 'threshold', 'niter', 'pblimit', 'pbmask', 'restfreq', 'nchan',
                           'width', 'start', 'chanchunks', 'deconvolver', 'weighting', 'robust', 'git_version', 'git_date', 'version']
    colnames_stats = ('min max std sum mean'.split() +
                      'min_K max_K std_K sum_K mean_K'.split() +
                      'lowmin lowmax lowstd lowmadstd lowsum lowmean'.split() +
                      ['mod' + x for x in 'min max std sum mean'.split()] + ['epsilon'])

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
        tbl.write(tbldir / 'cube_stats.ecsv', overwrite=True)
        tbl.write(tbldir / 'cube_stats.ipac', format='ascii.ipac', overwrite=True)
        tbl.write(tbldir / 'cube_stats.html', format='ascii.html', overwrite=True)
        tbl.write(tbldir / 'cube_stats.tex', overwrite=True)
        tbl.write(tbldir / 'cube_stats.js.html', format='jsviewer', overwrite=True)
        return tbl

    if os.getenv('START_FROM_CACHED') == 'False':
        start_from_cached = False  # TODO: make a parameter
    else:
        start_from_cached = True
        print(f"Starting from cached file {tbldir / 'cube_stats.ecsv'}")
    tbl = None
    if start_from_cached and os.path.exists(tbldir / 'cube_stats.ecsv'):
        tbl = Table.read(tbldir / 'cube_stats.ecsv')
        print(tbl)
        if len(tbl.colnames) != NCOLS:
            warnings.warn("Cached file is BAD!  Moving it.")
            shutil.move(tbldir / 'cube_stats.ecsv',
                        tbldir / f'cube_stats_{datetime.datetime.now().isoformat()}.ecsv')
            rows = []
        else:
            rows = [[tbl[cn].quantity[ii]
                     if tbl[cn].unit not in (None, u.dimensionless_unscaled)
                     else tbl[cn][ii] for cn in tbl.colnames]
                    for ii in range(len(tbl))]
    else:
        rows = []

    cache_stats_file = open(tbldir / "cube_stats.txt", 'w')

    mousmap = get_mous_to_sb_mapping('2021.1.00172.L')
    mousmap_ = {key.replace("/", "_").replace(":", "_"): val for key, val in mousmap.items()}

    for fullpath in glob.glob(f"{basepath}/sci*/group*/member*/"):
        mous = os.path.basename(fullpath.strip('/')).split(".")[-1]
        sbname = mousmap_[mous]
        field = sbname.split("_")[3]
        config = sbname.split("_")[5]
        if ' ' in config:
            # handle this case: 'Sgr_A_st_aj_03_7M Sgr_A_st_aj_03_7M_original'
            config = config.split(" ")[0]
        rerun = 'original' in sbname

        for suffix in (".image", ):  # ".contsub.image"):#, ".contsub.JvM.image.fits", ".JvM.image.fits"):
            globblob = f'{fullpath}/calibrated/working/*.iter1{suffix}'
            fns = glob.glob(globblob)
            globblob2 = f'{fullpath}/reclean/*.iter1{suffix}'
            fns += glob.glob(globblob2)

            for fn in fns:
                if 'mfs' in fn:
                    continue

                spw = int([x.lstrip('spw') for x in fn.split(".") if 'spw' in x][0])

                if tbl is not None:
                    row_matches = ((tbl['Field'] == field) &
                                   (tbl['Config'] == config) &
                                   (tbl['spw'] == spw) &
                                   (tbl['suffix'] == suffix))
                    if any(row_matches):
                        print(f"Skipping {globblob} as complete: {tbl[row_matches]}", flush=True)
                        continue

                modfn = fn.replace(".image", ".model")
                if os.path.exists(fn) and not os.path.exists(modfn):
                    log.error(f"File {fn} is missing its model {modfn}")
                    continue
                psffn = fn.replace(".image", ".psf")

                print(f"Beginning field {field} config {config} spw {spw} suffix {suffix}", flush=True)
                print(f"File: '{fn}'", flush=True)

                logtable = casaTable.read(f'{fn}/logtable')
                hist = logtable['MESSAGE']

                history = {x.split(":")[0]: ":".join(x.split(": ")[1:])
                           for x in hist if ':' in x and 'ICRS' not in x}
                history.update({x.split("=")[0]: x.split("=")[1].lstrip()
                                for x in hist if '=' in x})

                jvmimage = fn.replace(".image", ".JvM.image")
                if os.path.exists(jvmimage):
                    fn = jvmimage
                elif os.path.exists(jvmimage + ".fits"):
                    fn = jvmimage + ".fits"
                elif os.path.exists(fn):
                    pass
                elif os.path.exists(fn + ".fits"):
                    fn = fn + ".fits"

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

                if 'imsize' not in history:
                    history['imsize'] = str(cube.shape[1:])
                if 'cell' not in history:
                    history['cell'] = str([x.to(u.arcsec).to_string() for x in cube.wcs.celestial.proj_plane_pixel_scales()])
                if 'restfreq' not in history:
                    history['restfreq'] = float(cube.wcs.wcs.restfrq)
                if 'nchan' not in history:
                    history['nchan'] = int(cube.shape[0])

                with sched:
                    # mask to select the channels with little/less emission
                    meanspec = cube.mean(axis=(1, 2))
                    lowsignal = meanspec < np.nanpercentile(meanspec, 25)

                    print(f"Low-signal region selected {lowsignal.sum()} channels out of {lowsignal.size}."
                          f" ({lowsignal.sum() / lowsignal.size * 100:0.2f}) %")

                    assert lowsignal.sum() > 0
                    assert lowsignal.sum() < lowsignal.size

                    if False:
                        # this is an open to-do item: we need to create noise estimation regions
                        noiseregion = get_noise_region(field)
                        dt(f"Getting noise region {noiseregion}")
                        assert noiseregion is not None
                        noiseest_cube = cube.subcube_from_regions(regions.Regions.read(noiseregion))
                    noiseest_cube = cube

                    dt(cube)
                    dt(noiseest_cube)

                    minfreq = cube.spectral_axis.min()
                    maxfreq = cube.spectral_axis.max()
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

                #min = cube.min()
                #max = cube.max()
                ##mad = cube.mad_std()
                #std = cube.std()
                #sum = cube.sum()
                #mean = cube.mean()

                del stats
                del faintstats

                if os.path.exists(modfn):
                    modcube = SpectralCube.read(modfn, format='casa_image', target_chunksize=target_chunksize)
                elif os.path.exists(modfn + ".fits"):
                    modcube = SpectralCube.read(modfn + ".fits", format='fits', use_dask=True)
                modsched = modcube.use_dask_scheduler(scheduler=scheduler, num_workers=num_workers)

                dt(modcube)
                dt(f"Computing model cube statistics with scheduler {scheduler} and sched args {modcube._scheduler_kwargs}")
                with modsched:
                    modstats = modcube.statistics()
                dt("Done with model stats")
                modmin = modstats['min']
                modmax = modstats['max']
                modstd = modstats['sigma']
                modsum = modstats['sum']
                modmean = modstats['mean']

                del modcube
                del modstats

                if os.path.exists(psffn):
                    try:
                        (residual_peak, peakloc_as, frac, epsilon, firstnull, r_sidelobe, _) = get_psf_secondpeak(psffn, specslice=slice(cube.shape[0] // 2, cube.shape[0] // 2 + 1))
                    except Exception as ex:
                        print(f"get_psf_secondpeak failed with {ex} for file {psffn}")
                        continue

                del cube

                jtok_equiv = beam.jtok_equiv(u.Quantity(minfreq + maxfreq, u.Hz) / 2)

                row = ([field, config, spw, suffix, os.path.basename(fn), beam.major.to(u.arcsec), beam.minor.to(u.arcsec), beam.pa,
                        u.Quantity(restfreq, u.Hz), u.Quantity(minfreq, u.Hz), u.Quantity(maxfreq, u.Hz)] +
                       [history[key] if key in history else '' for key in colnames_fromheader] +
                       [min, max, std, sum, mean] +
                       list(map(lambda x: u.Quantity(x).to(u.K, jtok_equiv), [min, max, std, sum, mean])) +
                       [lowmin, lowmax, lowstd, lowmadstd, lowsum, lowmean] +
                       [modmin, modmax, modstd, modsum, modmean, epsilon])
                assert len(row) == len(colnames) == NCOLS
                rows.append(row)

                cache_stats_file.write(" ".join(map(str, row)) + "\n")
                cache_stats_file.flush()
                print(f'len(rows): {len(rows)}, len(colnames): {len(colnames)}')
                tbl = save_tbl(rows, colnames)

    cache_stats_file.close()

    print(tbl)

    os.chdir(cwd)

    if threads and nthreads > 1 and 'client' in locals():
        client.close()  # noqa
