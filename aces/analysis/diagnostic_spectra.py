import numpy as np
import json
import os
import radio_beam
import reproject
from astropy import constants, units as u, table, stats, coordinates, wcs, log
from astropy.io import fits
from spectral_cube import SpectralCube, wcs_utils, tests, Projection, OneDSpectrum
from astropy.nddata import Cutout2D
from aces.analysis.parse_contdotdat import parse_contdotdat
from aces.analysis import continuum_selection_diagnostic_plots
from aces import conf
import glob

import tempfile

import pylab as pl
pl.ioff()

overwrite=False

basepath = conf.basepath

def make_diagnostic_spectra(fn):
    basedir = os.path.dirname(fn)
    basename = os.path.splitext(os.path.basename(fn))[0]
    memdir = os.path.split(os.path.split(basedir)[0])[0]

    specdir = os.path.join(basedir, 'spectra')
    if not os.path.exists(specdir):
        os.mkdir(specdir)
        os.mkdir(os.path.join(specdir, 'pngs'))
    if os.path.exists(fn):
        cube = SpectralCube.read(fn,
                                    use_dask=True).with_spectral_unit(u.GHz)
    else:
        log.exception("File {0} does not exist".format(fn))
        return

    for operation in ('mean', 'max', 'median'):
        out_fn = f'{specdir}/{basename}.{operation}spec.fits'
        if overwrite or not os.path.exists(out_fn):
            spec = getattr(cube, operation)(axis=(1,2))
            spec.write(out_fn, overwrite=overwrite)

        spec_jy = OneDSpectrum.from_hdu(fits.open(out_fn)).with_spectral_unit(u.GHz)
        if cube.shape[0] != spec_jy.size:
            spec_jy = getattr(cube, operation)(axis=(1,2))
            spec_jy.write(out_fn, overwrite=True)

        try:
            jtok = cube.jtok_factors()
        except AttributeError:
            jtok = cube.beam.jtok(cube.spectral_axis)
        spec_K = spec_jy * jtok*u.K / (u.Jy/u.beam)

        for spec, unit in zip((spec_jy, spec_K), ("", "K")):
            fig_fn = f'{specdir}/pngs/{basename}.{operation}spec.png'
            pl.clf()
            spec.quicklook(filename=fig_fn, color='k',
                            linewidth=0.9, drawstyle='steps-mid')
            sel = np.zeros(spec.size, dtype='int')

            cdatfile = os.path.join(memdir, 'calibration/cont.dat')
            if os.path.exists(cdatfile):
                contfreqs = parse_contdotdat(cdatfile)

                for freqrange in contfreqs.split(";"):
                    low,high = freqrange.split("~")
                    high = u.Quantity(high)
                    low = u.Quantity(low, unit=high.unit)
                    sel += (spec.spectral_axis > low) & (spec.spectral_axis < high)
                    #print(f"{field}_{spw}: {low}-{high} count={sel.sum()}")

                usel = np.unique(sel)
                # 0 means 'not included in any windows', 1 means 'included in 1 window'
                # 2 or more means included in 2 or more.
                # The cases addressed here are:
                # {0,1}: some continuum, some not
                # {1}: all continuum
                # {0,1,2,...} or {1,2,...}: some or all continuum, at least one pixel twice or more
                if set(usel) in ({0, 1}, {1}):
                    sel = sel.astype('bool')

                    dat_to_plot = spec.value.copy()
                    dat_to_plot[~sel] = np.nan
                    pl.plot(spec.spectral_axis, dat_to_plot, linewidth=4,
                            zorder=-5, alpha=0.75, color='orange')
                elif len(usel) > 1:
                    dat_to_plot = np.empty(spec.value.shape)
                    dat_to_plot[:] = np.nan
                    # skip zero
                    for selval in usel[1:]:
                        dat_to_plot[sel == selval] = spec.value[sel == selval]
                    pl.plot(spec.spectral_axis, dat_to_plot, linewidth=4,
                            zorder=selval-10, alpha=0.75, color='orange')
                else:
                    log.error(f"No selected continuum for {basename}.{operation}: {sel.sum()} {usel}")
                    continue
                print(f"{basename}: {sel.sum()} {usel}")
            pl.title(f"{basename} {operation}")
            pl.savefig(fig_fn, bbox_inches='tight')


def main():
    # need to be in main block for dask to work
    #from dask.distributed import Client
    #if os.getenv('SLURM_MEM_PER_NODE'):
    #    memlim_total = int(os.getenv('SLURM_MEM_PER_NODE')) / 1024 # GB
    #    ntasks = int(os.getenv('SLURM_NTASKS'))
    #    memlim = memlim_total / ntasks
    #else:
    #    memlim = 1
    #    ntasks = 8
    #client = Client(memory_limit=f'{memlim}GB', n_workers=ntasks)
    #nworkers = len(client.scheduler_info()['workers'])
    #print(f"Client schedular info: {client.scheduler_info()['services']}")
    #print(f"Number of workers: {nworkers}")
    #print(f"Client schedular info: {client.scheduler_info()}")
    #print(f"Client vers: {client.get_versions(check=True)}")
    if os.getenv('ENVIRONMENT') == 'BATCH':
        pass
    else:
        from dask.diagnostics import ProgressBar
        pbar = ProgressBar()
        pbar.register()

    gpath = 'data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9'
    
    files = sorted(glob.glob(f'{basepath}/{gpath}/member*/calibrated/working/*.statcont.contsub.fits'))

    if os.getenv('SLURM_ARRAY_TASK_ID') is not None:
        slurm_array_task_id = int(os.getenv('SLURM_ARRAY_TASK_ID'))
    else:
        slurm_array_task_id = None

    for ii, fn in enumerate(files):
        # either if the task ID is specified and matches this one, or if it's unspecified
        if slurm_array_task_id in (ii, None):
            if slurm_array_task_id is not None:
                print(ii, fn)
            make_diagnostic_spectra(fn)

            memberid = "_".join(fn.split("/")[-4].split("_")[-2:])
            continuum_selection_diagnostic_plots.make_plot(memberid)


def get_file_numbers():
    """
    For slurm jobs, just run through all the files that we're maybe going to make diagnostic spectra for and check which ones need it
    """

    redo = bool(os.getenv('REDO'))

    filenames = sorted(glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_*/calibrated/working/*statcont.contsub.fits'))

    numlist = []
    for ii, fn in enumerate(filenames):

        basename = os.path.splitext(os.path.basename(fn))[0]
        basedir = os.path.dirname(fn)
        specdir = os.path.join(basedir, 'spectra')
        for operation in ('max', 'mean', 'median'):
            out_fn = f'{specdir}/{basename}.{operation}spec.fits'

            if not os.path.exists(out_fn) or redo:
                numlist.append(ii)
    
    return numlist



if __name__ == "__main__":
    main()