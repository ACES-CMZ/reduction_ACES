"""
Continuum-included and positive-only moment-0 maps of a giant line cube (CS 2-1 by default).

Continuum subtraction changes the line integral wherever there is continuum
behind the gas.  For gas with excitation temperature T_ex and optical depth
tau in front of continuum T_C,

    T_obs   = T_C exp(-tau) + T_ex (1 - exp(-tau))
    Delta T = T_obs - T_C = (T_ex - T_C) (1 - exp(-tau))

so the moment 0 of the continuum-subtracted cube is the line emission minus
the continuum absorbed by the gas, (1 - T_C / T_ex) times the line emission.
Toward bright continuum (T_C ~ T_ex) the continuum-subtracted moment 0
underestimates the line, and where T_C > T_ex emitting gas appears as a
negative.  A compact continuum source beside (not behind) the gas in the same
beam gives the same spectrum; see the figure in
https://github.com/ACES-CMZ/reduction_ACES/pull/476 .

The giant cubes are built from continuum-subtracted cubes.  This script adds
the reimaged continuum back, convolved to the cube's beam, and writes moment-0
maps within the line's velocity window and within the persisted
``signal_mask_pruned`` (``hlsig_dilated``) mask:

- ``withcont``          : moment 0 of (cube + continuum), i.e. of T_obs
- ``positive``          : moment 0 of the cube with negative voxels excluded
- ``withcont_positive`` : moment 0 of (cube + continuum) with negative voxels excluded

The positive-only maps exclude absorption, real or produced by continuum
subtraction.  The convolved continuum map is also written, so that
T_C / T_ex corrections can be made from the same continuum.

Usage::

    MOLNAME=CS21 python giantcube_contincluded_moments.py

``CONTINUUM_FILE`` overrides the continuum mosaic (a filename in
``{basepath}/mosaics/continuum/``).
"""
import os
import time

import numpy as np
import pylab as pl
import radio_beam
from astropy import units as u
from astropy.io import fits
from dask.config import set as dask_set
from dask.diagnostics import ProgressBar
from spectral_cube import SpectralCube, Projection, BooleanArrayMask

from aces import conf
from aces.imaging.make_mosaic import makepng
from aces.analysis.giantcube_cuts import velocity_mask, update_mom0

basepath = conf.basepath
cubepath = f'{basepath}/mosaics/cubes/'
mompath = f'{basepath}/mosaics/cubes/moments/'
contpath = f'{basepath}/mosaics/continuum/'

# per-line configuration: the reimaged continuum mosaic that covers the line's
# spectral window, and the velocity window used by do_all_stats in giantcube_cuts
LINES = {
    # CS 2-1 (97.98 GHz) is in spw33
    'CS21': dict(continuum='12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits',
                 vmask_kwargs=dict(vmax_pos=300 * u.km / u.s, vmax_neg=200 * u.km / u.s,
                                   vmin_pos=-200 * u.km / u.s, vmin_neg=-300 * u.km / u.s)),
}


def continuum_at_cube_beam(cube, contfn):
    """Read the continuum mosaic and convolve it to the cube's beam.

    The continuum mosaics share the giant cubes' celestial grid.  If the
    continuum beam cannot be deconvolved from the cube beam (the reimaged
    continuum beam is circular; the cube beam is elliptical with a slightly
    smaller minor axis), the continuum is convolved to the smallest beam
    enclosing both.
    """
    cont = Projection.from_hdu(fits.open(contfn)[0])
    assert cont.shape == cube.shape[1:], f"{contfn} is not on the cube's grid"
    commonbeam = radio_beam.Beams(beams=[cube.beam, cont.beam]).common_beam()
    print(f"Continuum beam {cont.beam}; cube beam {cube.beam}; convolving continuum to {commonbeam}",
          flush=True)
    # convolve_to rescales Jy/beam images by the beam-area ratio
    return cont.convolve_to(commonbeam).to(cube.unit)


def contincluded_moments(cube, molname, cont, mask, howargs={}):
    """Write withcont / positive / withcont_positive moment-0 maps of ``cube``
    within the velocity window and within ``mask``."""
    cfg = LINES[molname]
    t0 = time.time()

    cube = cube.with_mask(velocity_mask(cube, **cfg['vmask_kwargs']))
    fwidth = cube.spectral_axis.max() - cube.spectral_axis.min()
    nocont = ~np.isfinite(cont.value)

    for maskname, mcube in (('', cube), ('masked_hlsig_dilated_', cube.with_mask(mask))):
        # plain quantity: the continuum beam can differ slightly from the cube beam
        withcont = mcube + u.Quantity(cont)
        products = {'withcont': withcont,
                    'positive': mcube.with_mask(mcube > 0 * mcube.unit),
                    'withcont_positive': withcont.with_mask(withcont > 0 * withcont.unit)}
        for pname, pcube in products.items():
            print(f"{maskname}{pname} mom0.  dt={time.time() - t0:0.1f}", flush=True)
            mom0 = pcube.moment0(axis=0, **howargs).to(u.Jy / u.beam * u.km / u.s)
            if 'withcont' in pname:
                mom0[nocont] = np.nan
            mom0path = f"{mompath}/{molname}_CubeMosaic_{maskname}{pname}_mom0.fits"
            mom0.write(mom0path, overwrite=True)
            update_mom0(mom0path=mom0path, fwidth=fwidth, cube=cube)
            with fits.open(mom0path, mode='update') as hdul:
                if 'withcont' in pname:
                    hdul[0].header['HISTORY'] = f'continuum added from {cfg["continuum"]}'
                if 'positive' in pname:
                    hdul[0].header['HISTORY'] = 'negative voxels excluded from the sum'
            makepng(data=mom0.value, wcs=mom0.wcs, imfn=mom0path.replace('.fits', '.png'),
                    stretch='asinh', vmin=-0.1, max_percent=99.5)
            print(f"Wrote {mom0path}.  dt={time.time() - t0:0.1f}", flush=True)


def main():
    molname = os.getenv('MOLNAME', 'CS21')
    contfile = os.getenv('CONTINUUM_FILE', LINES[molname]['continuum'])
    LINES[molname]['continuum'] = contfile

    nworkers = int(os.getenv('SLURM_NTASKS', os.getenv('SLURM_CPUS_PER_TASK', 1)))
    dask_set({"array.chunk-size": "1GB"})

    cube = SpectralCube.read(f'{cubepath}/{molname}_CubeMosaic.fits', use_dask=True)
    cube.use_dask_scheduler('threads', num_workers=nworkers)
    ProgressBar().register()
    print(cube, flush=True)

    cont = continuum_at_cube_beam(cube, f'{contpath}/{contfile}')
    contoutfn = f'{mompath}/{molname}_CubeMosaic_continuum.fits'
    cont.write(contoutfn, overwrite=True)
    with fits.open(contoutfn, mode='update') as hdul:
        hdul[0].header['HISTORY'] = f'{contfile} convolved to the {molname} cube beam'
    makepng(data=cont.value, wcs=cont.wcs, imfn=contoutfn.replace('.fits', '.png'),
            stretch='asinh', min_percent=1, max_percent=99.5, cmap=pl.cm.gray_r)

    # persisted 3D signal mask from giantcube_cuts.do_all_stats (values 0/1).  Its
    # header has only celestial WCS keys, so it is given the cube's WCS before reading.
    maskhdu = fits.open(f'{mompath}/{molname}_CubeMosaic_signal_mask_pruned.fits', memmap=True)[0]
    maskhdu.header.update(cube.wcs.to_header())
    maskcube = SpectralCube.read(maskhdu, use_dask=True)
    assert maskcube.shape == cube.shape
    mask = BooleanArrayMask((maskcube._data > 0.5).rechunk(cube._data.chunks), cube.wcs)

    contincluded_moments(cube, molname, cont, mask)


if __name__ == "__main__":
    main()
