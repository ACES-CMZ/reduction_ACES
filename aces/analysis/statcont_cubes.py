"""
Ported from ALMA-IMF

It is intended to be run in `imaging_results/` and will produce statcont
contfiles

It is partly a performance test - for the bigger cubes, there were sometimes
memory problems

The noise estimation region and num_workers are both hard-coded and should be
customized.

We should eventually allow multi-cube combination using full statcont abilities
"""
import time
import warnings
from astropy.table import Table
from spectral_cube import SpectralCube
import spectral_cube
from spectral_cube.utils import NoBeamError
from astropy.io import fits
import dask

from statcont.cont_finding import c_sigmaclip_scube

import glob

import tempfile

import os
import sys

from aces import conf
from aces.imaging.make_mosaic import makepng

basepath = conf.basepath

# for zarr storage
if not os.getenv('SLURM_TMPDIR'):
    os.environ['TMPDIR'] = '/blue/adamginsburg/adamginsburg/tmp'


def main():
    # need to be in main block for dask to work
    #from dask.distributed import Client
    #if os.getenv('SLURM_MEM_PER_NODE'):
    #    memlim_total = int(os.getenv('SLURM_MEM_PER_NODE')) / 1024 # GB
    #    ntasks = int(os.getenv('SLURM_NTASKS'))
    #    memlim = memlim_total / ntasks
    #    print(f"Memory limit is {memlim} GB")
    #else:
    #    memlim = 1
    #    ntasks = 8
    #client = Client(memory_limit=f'{memlim}GB', n_workers=ntasks)
    #nworkers = len(client.scheduler_info()['workers'])
    #print(f"Client scheduler info: {client.scheduler_info()['services']}")
    #print(f"Number of workers: {nworkers}  (should be equal to ntasks={ntasks})")
    #print(f"Client scheduler info: {client.scheduler_info()}")
    #print(f"Client vers: {client.get_versions(check=True)}")
    if os.getenv('ENVIRONMENT') == 'BATCH':
        pass
    else:
        from dask.diagnostics import ProgressBar
        pbar = ProgressBar()
        pbar.register()

    nthreads = os.getenv('SLURM_NTASKS')
    if nthreads is not None:
        nthreads = int(nthreads)
        dask.config.set(scheduler='threads')
    else:
        dask.config.set(scheduler='synchronous')

    scheduler = dask.config.get('scheduler')
    print(f"Using {nthreads} threads with the {scheduler} scheduler")

    print(f"tempdir is {tempfile.gettempdir()}")

    # kinda meant to be hacked
    redo = False

    tbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/cube_stats.ecsv')

    def get_size(start_path='.'):
        total_size = 0
        for dirpath, dirnames, filenames in os.walk(start_path):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                # skip if it is symbolic link
                if not os.path.islink(fp):
                    total_size += os.path.getsize(fp)

        return total_size

    # simpler approach
    #sizes = {fn: get_size(fn) for fn in glob.glob(f"{basepath}/*_12M_spw[0-9].image")}
    #filenames = [f'{basepath}/{fn}' for fn in tbl['filename']] + list(glob.glob(f"{basepath}/*_12M_spw[0-9].image")) + list(glob.glob(f"{basepath}/*_12M_sio.image"))
    filenames = glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_*/calibrated/working/*cube*.image.pbcor.fits')

    sizes = {ii: get_size(fn)
             for ii, fn in enumerate(filenames)
             }

    target_chunk_size = int(1e8)

    if os.getenv('SLURM_ARRAY_TASK_ID') is not None:
        slurm_array_task_id = int(os.getenv('SLURM_ARRAY_TASK_ID'))
    else:
        slurm_array_task_id = None

    for ii in sorted(sizes, key=lambda x: sizes[x]):

        if slurm_array_task_id is not None:
            if ii != slurm_array_task_id:
                continue

        fn = filenames[ii]

        #outfn = fn+'.statcont.cont.fits'
        outfn = fn.replace(".image.pbcor.fits", ".image.pbcor.statcont.cont.fits")
        fileformat = 'fits'
        assert outfn.count('.fits') == 1

        if not os.path.exists(outfn) or redo:
            t0 = time.time()

            # touch the file to allow parallel runs
            with open(outfn, 'w') as fh:
                fh.write("")

            print(f"{fn}->{outfn}, size={sizes[ii] / 1024**3} GB", flush=True)

            print(f"Target chunk size is {target_chunk_size}", flush=True)
            cube = SpectralCube.read(fn, target_chunk_size=target_chunk_size,
                                     format=fileformat, use_dask=True)
            if 'JvM' not in fn:
                print(f"Minimizing {cube}", flush=True)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    cube = cube.minimal_subcube()
            print(cube, flush=True)
            sys.stdout.flush()
            sys.stderr.flush()

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                print(f"Doing statcont with {nthreads} threads")
                with cube.use_dask_scheduler('threads', num_workers=nthreads):
                    print("Calculating noise", flush=True)
                    if ii < len(tbl):
                        noise = tbl['std'].quantity[ii]
                    else:
                        noise = cube.std()

                    print("Sigma clipping", flush=True)
                    result = c_sigmaclip_scube(cube, noise,
                                               verbose=True,
                                               save_to_tmp_dir=True)
                    print("Running the compute step", flush=True)
                    data_to_write = result[1].compute()
                    cont = data_to_write.value

                    print(f"Writing to FITS {outfn}", flush=True)
                    fits.PrimaryHDU(data=cont,
                                    header=cube[0].header).writeto(outfn,
                                                                   overwrite=True)
            print(f"{fn} -> {outfn} in {time.time() - t0}s", flush=True)
        else:
            try:
                cont = fits.getdata(outfn)
            except Exception as ex:
                shutil.move(outfn, outfn.replace(".fits", ".bad.fits"))
                print(f"File {outfn} exists but could not be opened; renaming it.  Try again.", flush=True)
                print(ex)
                continue
            print(f"{fn} is done, loaded {outfn}", flush=True)

        if os.path.exists(outfn):
            if os.path.getsize(outfn) == 0:
                print(f"{outfn} had size {os.path.getsize(outfn)}", flush=True)
                os.remove(outfn)

        outcube = contsubfn = fn.replace(".image.pbcor.fits", ".image.pbcor.statcont.contsub.fits")
        if os.path.exists(contsubfn):
            try:
                print(f"Checking {contsubfn} for beam")
                beam = SpectralCube.read(contsubfn).beam
            except NoBeamError:
                try:
                    beam = SpectralCube.read(fn).beam
                except NoBeamError:
                    raise ValueError(f"Neither {fn} nor {contsubfn} have a beam")
                except AttributeError:
                    beams = SpectralCube.read(fn).beams
                    raise AttributeError(f"{fn} is a multi-beam cube")
                redo = True
            except AttributeError:
                beams = SpectralCube.read(contsubfn).beams
                redo = True

        if fn.endswith('.fits'):
            assert outcube.count('.fits') == 1
            if (not os.path.exists(outcube)) or redo:
                print(f"Writing contsub cube to {outcube}", flush=True)
                cube = SpectralCube.read(fn,
                                         target_chunk_size=target_chunk_size,
                                         use_dask=True, format=fileformat)
                cube.allow_huge_operations = True
                if cube.shape[1] != cont.shape[0] or cube.shape[2] != cont.shape[1]:
                    print(f"Minimizing {cube}", flush=True)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        cube = cube.minimal_subcube()
                    cube.allow_huge_operations = True

                scube = cube - cont * cube.unit
                scube.write(outcube, overwrite=True)
            else:
                print(f"Found existing cube {outcube} and redo=False")
        else:
            print("Operated on a non-FITS file: no contsub cube was created")
        sys.stdout.flush()
        sys.stderr.flush()


if __name__ == "__main__":
    main()