"""
Continuum-included and positive-only moment-0 maps of a giant line cube (CS 2-1 by default).

The giant cubes are mosaics of feathered (TP + 7m + 12m) cubes built from
statcont continuum-subtracted 7m and 12m cubes.  statcont over-subtracts the
continuum where the line emission is bright and broad, producing artificial
absorption and biased integrated intensities.

This script:

1. Rebuilds the continuum that statcont removed: for each field, the 12m and 7m
   ``*.statcont.cont.fits`` images are feathered together (TP line cubes are
   baseline-subtracted, so they carry no continuum), convolved to the giant
   cube's common beam, and mosaicked onto the giant-cube grid with the same 12m
   weights used to build the cube.  Adding this map to the cube restores the
   line + continuum cube that went into statcont.

2. Makes moment-0 maps within the line's velocity window (``velocity_mask``)
   and within the persisted ``signal_mask_pruned`` (the ``hlsig_dilated`` mask):

   - ``withcont``          : sum of (cube + continuum)
   - ``positive``          : sum of max(cube, 0)
   - ``withcont_positive`` : sum of max(cube + continuum, 0)

   The positive-only maps exclude absorption (real or artificial).

Usage::

    MOLNAME=CS21 python giantcube_contincluded_moments.py

Set ``SKIP_CONTINUUM=True`` to reuse an existing continuum mosaic.
"""
import glob
import os
import time

import numpy as np
import pandas as pd
import radio_beam
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from spectral_cube import SpectralCube, Projection

from aces import conf
from aces.imaging.make_mosaic import make_mosaic, makepng, _weight_is_flat
from aces.imaging.mosaic_12m import get_weightfile
from aces.analysis.giantcube_cuts import velocity_mask

basepath = conf.basepath
cubepath = f'{basepath}/mosaics/cubes/'
mompath = f'{basepath}/mosaics/cubes/moments/'
rawpath = f'{basepath}/rawdata/2021.1.00172.L/'

# per-line configuration: 12m spw, 7m spw, feathered-cube glob, velocity window
LINES = {
    'CS21': dict(spw12=33, spw7=24, featherdir='SPW33', featherspw='SPW_33',
                 vmask_kwargs=dict(vmax_pos=300 * u.km / u.s, vmax_neg=200 * u.km / u.s,
                                   vmin_pos=-200 * u.km / u.s, vmin_neg=-300 * u.km / u.s)),
}

# preference order for which statcont continuum to use when there are several
# imaging versions; matches the preference in get_weightfile
VERSION_PREFERENCE = ('iter1.reclean', 'iter1', 'manual', 'manual.reclean', 'iter2.reclean')


def find_statcont_cont(mousid, spw):
    """Find the statcont continuum image for one MOUS and spw."""
    candidates = glob.glob(f"{rawpath}/s*/g*/member.uid___A001_{mousid}/calibrated/working/"
                           f"*.Sgr_A_star_sci.spw{spw}.cube.I.*.image.pbcor.statcont.cont.fits")
    for version in VERSION_PREFERENCE:
        matches = [fn for fn in candidates
                   if f'.cube.I.{version}.image.pbcor.statcont.cont.fits' in fn]
        if matches:
            # s38 is the later pipeline run; prefer it over s12 duplicates
            return sorted(matches, key=lambda x: '.s38_' not in x)[0]
    return None


def find_fallback_cont(mousid, spw12):
    """For fields without a statcont continuum, use a directly imaged continuum
    covering the line's spw."""
    candidates = sorted(glob.glob(f"{rawpath}/s*/g*/member.uid___A001_{mousid}/calibrated/working/"
                                  f"*.Sgr_A_star_sci.spw{spw12}_*.cont.I.*image.tt0.pbcor.fits"))
    return candidates[0] if candidates else None


def read_2d(fn):
    hdu = fits.open(fn)[0]
    data = hdu.data.squeeze()
    header = WCS(hdu.header).celestial.to_header()
    beam = radio_beam.Beam.from_fits_header(hdu.header)
    header.update(beam.to_header_keywords())
    header['BUNIT'] = hdu.header['BUNIT']
    return fits.PrimaryHDU(data=data, header=header)


def feather_cont(hdu12, hdu7):
    """Feather the 7m continuum into the 12m continuum, keeping the 12m footprint."""
    import uvcombine
    nan12 = ~np.isfinite(hdu12.data)
    hi = fits.PrimaryHDU(data=np.nan_to_num(hdu12.data), header=hdu12.header)
    lo = fits.PrimaryHDU(data=np.nan_to_num(hdu7.data), header=hdu7.header)
    feathered = uvcombine.feather_simple(hi, lo, return_hdu=True)
    data = np.real(feathered.data)
    data[nan12] = np.nan
    return fits.PrimaryHDU(data=data, header=hdu12.header)


def weight_plane(fn):
    """Read a 2D plane of a 12m weight file (mid-channel for a weight cube)."""
    if _weight_is_flat(fn):
        if fn.endswith('.fits'):
            prj = Projection.from_hdu(fits.open(fn)[0])
        else:
            prj = SpectralCube.read(fn, format='casa_image')[0]
    else:
        wcube = SpectralCube.read(fn, format='fits' if fn.endswith('fits') else 'casa_image',
                                  use_dask=True)
        prj = wcube[wcube.shape[0] // 2]
    hdu = prj.hdu
    hdu.header = WCS(hdu.header).celestial.to_header()
    return hdu


def make_continuum_mosaic(molname, outfn):
    cfg = LINES[molname]
    sb = pd.read_csv(f'{basepath}/reduction_ACES/aces/data/tables/aces_SB_uids.csv')
    cubehead = fits.getheader(f'{cubepath}/{molname}_CubeMosaic.fits')
    target_header = WCS(cubehead).celestial.to_header()
    target_header['NAXIS'] = 2
    target_header['NAXIS1'] = cubehead['NAXIS1']
    target_header['NAXIS2'] = cubehead['NAXIS2']
    commonbeam = radio_beam.Beam.from_fits_header(cubehead)

    featherfiles = glob.glob(f'{basepath}/upload/Feather_12m_7m_TP/{cfg["featherdir"]}/cubes/'
                             f'Sgr_A_st_*.TP_7M_12M_feather_all.{cfg["featherspw"]}.image.statcont.contsub.fits')
    fields_in_cube = {fn.split('Sgr_A_st_')[-1].split('.')[0]: fn for fn in featherfiles}

    hdus, weights, notes = [], [], []
    for _, row in sb.iterrows():
        field = row['Obs ID']
        if field not in fields_in_cube:
            notes.append(f"{field}: not in giant cube; skipped")
            continue
        fn12 = find_statcont_cont(row['12m MOUS ID'], cfg['spw12'])
        fn7 = find_statcont_cont(row['7m MOUS ID'], cfg['spw7'])
        if fn12 is None:
            fn12 = find_fallback_cont(row['12m MOUS ID'], cfg['spw12'])
            if fn12 is None:
                notes.append(f"{field}: NO 12m continuum found; field omitted")
                continue
            notes.append(f"{field}: no 12m statcont continuum; using fallback {os.path.basename(fn12)}")
        hdu12 = read_2d(fn12)
        if fn7 is not None:
            hdu = feather_cont(hdu12, read_2d(fn7))
            notes.append(f"{field}: 12m={os.path.basename(fn12)} 7m={os.path.basename(fn7)}")
        else:
            hdu = hdu12
            notes.append(f"{field}: 12m={os.path.basename(fn12)} 7m=NONE (12m only)")

        wfn = get_weightfile(fields_in_cube[field], spw=cfg['spw12'])
        whdu = weight_plane(wfn)
        # put the weight on the continuum image's grid
        wdata, _ = reproject_interp(whdu, hdu.header, shape_out=hdu.data.shape)
        wdata[~np.isfinite(wdata) | ~np.isfinite(hdu.data)] = 0
        # same footprint cut as the giant cube (min_weight_fraction=0.05)
        wdata[wdata < 0.05 * np.nanmax(wdata)] = 0

        # convolve to the giant cube's common beam; in Jy/beam the values scale
        # with the beam area (same convention as make_mosaic)
        prj = Projection.from_hdu(hdu)
        beam = prj.beam
        cbhdu = prj.convolve_to(commonbeam).hdu
        cbhdu.data = cbhdu.data * (commonbeam.sr / beam.sr).decompose().value

        hdus.append(fits.HDUList([cbhdu]))
        weights.append(fits.PrimaryHDU(data=wdata, header=hdu.header))

    print("\n".join(notes), flush=True)

    # make_mosaic writes to {basepath}/mosaics/{folder}/{array}_{name}_mosaic.fits
    name = f'{molname}_statcont_continuum'
    make_mosaic(hdus, name=name, weights=weights, target_header=target_header,
                cbar_unit='Jy/beam', array='12m7m', folder='continuum',
                basepath=basepath, doplots=False)
    made = f'{basepath}/mosaics/continuum/12m7m_{name}_mosaic.fits'
    with fits.open(made, mode='update') as hdul:
        hdul[0].header.update(commonbeam.to_header_keywords())
        hdul[0].header['HISTORY'] = 'statcont continuum (12m feathered with 7m) per field,'
        hdul[0].header['HISTORY'] = 'convolved to giant-cube common beam, weighted by 12m weights'
        for line in notes:
            hdul[0].header['HISTORY'] = line[:70]
    if made != outfn:
        os.replace(made, outfn)
    return outfn


def contincluded_moments(molname, contfn, nrows=128):
    cfg = LINES[molname]
    cubefn = f'{cubepath}/{molname}_CubeMosaic.fits'
    maskfn = f'{mompath}/{molname}_CubeMosaic_signal_mask_pruned.fits'

    cube = SpectralCube.read(cubefn)
    vmask = velocity_mask(cube, **cfg['vmask_kwargs'])  # (nchan, 1, nx)
    dv = np.abs(cube.spectral_axis[1] - cube.spectral_axis[0]).to(u.km / u.s).value
    fwidth = cube.spectral_axis.max() - cube.spectral_axis.min()
    del cube

    cont = fits.getdata(contfn)
    # the giant cube is ~127 GB; memmap it and loop over strips of rows
    cubedata = fits.open(cubefn, memmap=True, ignore_missing_end=True)[0].data
    maskdata = fits.open(maskfn, memmap=True)[0].data
    nchan, ny, nx = cubedata.shape
    assert maskdata.shape == cubedata.shape
    assert cont.shape == (ny, nx)

    products = ['withcont', 'positive', 'withcont_positive']
    out = {(m, p): np.full((ny, nx), np.nan, dtype='float32')
           for m in ('vmask', 'hlsig') for p in products}

    t0 = time.time()
    for y0 in range(0, ny, nrows):
        y1 = min(y0 + nrows, ny)
        data = np.asarray(cubedata[:, y0:y1, :], dtype='float64')
        cstrip = cont[None, y0:y1, :]
        finite = np.isfinite(data)
        anyfinite = finite.any(axis=0)
        withcont = data + cstrip
        planes = {'withcont': withcont,
                  'positive': np.where(data > 0, data, 0),
                  'withcont_positive': np.where(withcont > 0, withcont, 0)}
        masks = {'vmask': vmask & finite,
                 'hlsig': vmask & finite & (np.asarray(maskdata[:, y0:y1, :]) > 0)}
        for mname, msk in masks.items():
            for pname, plane in planes.items():
                mom0 = np.where(msk, plane, 0).sum(axis=0) * dv
                # no-coverage pixels stay NaN; continuum-included maps need continuum coverage
                bad = ~anyfinite
                if 'withcont' in pname:
                    bad = bad | ~np.isfinite(cstrip[0])
                mom0[bad] = np.nan
                out[(mname, pname)][y0:y1, :] = mom0
        print(f"rows {y0}-{y1} of {ny} done. dt={time.time() - t0:0.1f}", flush=True)

    header = fits.getheader(f'{mompath}/{molname}_CubeMosaic_mom0.fits')
    header['CDELT3'] = (fwidth.value, 'Bandwidth integrated to make moment 0 map')
    outnames = {}
    for (mname, pname), arr in out.items():
        prefix = 'masked_hlsig_dilated_' if mname == 'hlsig' else ''
        outfn = f'{mompath}/{molname}_CubeMosaic_{prefix}{pname}_mom0.fits'
        hdr = header.copy()
        hdr['HISTORY'] = f'{pname} moment 0 ({mname} mask)'
        if 'withcont' in pname:
            hdr['HISTORY'] = f'continuum added from {os.path.basename(contfn)}'
        if 'positive' in pname:
            hdr['HISTORY'] = 'negative voxels set to zero before summing'
        fits.PrimaryHDU(data=arr, header=hdr).writeto(outfn, overwrite=True)
        makepng(data=arr, wcs=WCS(hdr).celestial, imfn=outfn.replace('.fits', '.png'),
                stretch='asinh', vmin=-0.1, max_percent=99.5)
        outnames[(mname, pname)] = outfn
        print(f"Wrote {outfn}", flush=True)
    return outnames


def main():
    molname = os.getenv('MOLNAME', 'CS21')
    contfn = f'{mompath}/{molname}_CubeMosaic_statcont_continuum.fits'
    skip = os.getenv('SKIP_CONTINUUM', 'False').lower() == 'true'
    if not (skip and os.path.exists(contfn)):
        make_continuum_mosaic(molname, contfn)
    contincluded_moments(molname, contfn)


if __name__ == "__main__":
    main()
