"""
Full-band, spatially-downsampled per-spw cubes via the zarr reproject_and_coadd
combiner (the "new approach").

For each spectral window (25, 27, 33, 35) this builds ONE cube covering the
*entire* spw bandwidth (all native channels), reprojected onto a coarse **5
arcsec** pixel grid and smoothed to a **10 arcsec** common beam.  These are
low-resolution full-band overview cubes -- cheap to make and to browse -- as
opposed to the per-line giant cubes.

Why the zarr 3D combiner rather than the per-channel `make_giant_mosaic_cube`:
the single `reproject.mosaicking.reproject_and_coadd(return_type='zarr')` call
does one coherent 3D reprojection, so it avoids the per-channel two_closest /
weight-slab machinery (and its data-vs-weight grid-offset NaN-striping bug).

Weights: the per-field weight cubes are essentially constant along frequency, so
we take each field's mid-channel weight plane, reproject it onto that field's
celestial grid, and broadcast it across the field's channels.  This gives a
correctly-located per-field weight without a full 3D weight reprojection and
without the shape-mismatch hazard of the per-channel path.

Resolution handling: fields are coadded at 5 arcsec sampling keeping their
native (~2.5 arcsec) beams, then the assembled 5 arcsec mosaic is convolved to
the 10 arcsec target beam (cheap at coarse sampling; 10 arcsec dominates the
per-field beam spread).  This avoids an expensive native-resolution convolution.

NB: uses the reproject `return_type='zarr'` feature (astropy/reproject#621) and
must be run in an environment where that reproject and a compatible
spectral-cube are installed (see aces/imaging/README_zarr_combine.md).  Set up
here but not yet run end-to-end -- treat as review-ready scaffolding.
"""
import glob
import os
import shutil

import numpy as np
import radio_beam
from astropy import units as u
from astropy import wcs
from astropy.io import fits

from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from spectral_cube import SpectralCube

from aces.imaging.mosaic_12m import get_weightfile
from aces import conf

basepath = conf.basepath

FEATHER = f'{basepath}/upload/Feather_12m_7m_TP'
NATIVE_HEADER = f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr'

# All ACES science spws we build full-band overview cubes for.
SPWS = (25, 27, 33, 35)


def spw_filelist(spw):
    return sorted(glob.glob(f'{FEATHER}/SPW{spw}/cubes/'
                            f'Sgr_A_st_*.TP_7M_12M_feather_all.SPW_{spw}.image.statcont.contsub.fits'))


def build_target_header(spw, ref_cube_file, pixscale=5 * u.arcsec,
                        native_header=NATIVE_HEADER):
    """3D target header: the native celestial footprint downsampled to
    ``pixscale`` pixels, plus the spw's native frequency axis taken from
    ``ref_cube_file``."""
    nat = fits.Header.fromtextfile(native_header)
    native_pix = abs(nat['CDELT2']) * u.deg
    factor = float((pixscale.to(u.deg) / native_pix).decompose())

    hdr = nat.copy()
    hdr['NAXIS'] = 3
    hdr['WCSAXES'] = 3
    # coarsen the two celestial axes by ``factor`` about the reference pixel
    for ax in (1, 2):
        n = nat[f'NAXIS{ax}']
        hdr[f'NAXIS{ax}'] = int(np.ceil(n / factor))
        hdr[f'CDELT{ax}'] = nat[f'CDELT{ax}'] * factor
        hdr[f'CRPIX{ax}'] = (nat[f'CRPIX{ax}'] - 0.5) / factor + 0.5

    # spectral axis: copy the spw's native FREQ grid
    refh = fits.getheader(ref_cube_file)
    refw = wcs.WCS(refh).spectral
    n3 = refh['NAXIS3']
    f0 = refw.pixel_to_world(0).to(u.Hz).value
    f1 = refw.pixel_to_world(n3 - 1).to(u.Hz).value
    hdr['NAXIS3'] = n3
    hdr['CTYPE3'] = 'FREQ'
    hdr['CUNIT3'] = 'Hz'
    hdr['CRPIX3'] = 1
    hdr['CRVAL3'] = f0
    hdr['CDELT3'] = (f1 - f0) / (n3 - 1)
    hdr['SPECSYS'] = 'LSRK'
    return hdr


def _weight_plane_on_grid(weightfile, refcube):
    """Reproject a field's mid-channel weight plane onto ``refcube``'s celestial
    grid and return it as a 2D array (weights are ~constant along frequency, so a
    single plane suffices; it is broadcast across channels by the caller)."""
    if weightfile.endswith('.fits'):
        whdu = fits.open(weightfile)[0]
        wdata = np.squeeze(whdu.data)
        if wdata.ndim == 3:
            wdata = wdata[wdata.shape[0] // 2]
        wwcs = wcs.WCS(whdu.header).celestial
    else:
        wc = SpectralCube.read(weightfile, format='casa_image', use_dask=True)
        wdata = np.squeeze(np.asarray(wc[wc.shape[0] // 2] if wc.ndim == 3 else wc))
        wwcs = wc.wcs.celestial
    plane, _ = reproject_interp((np.nan_to_num(wdata), wwcs),
                                refcube.wcs.celestial, shape_out=refcube.shape[1:])
    return np.nan_to_num(plane)


def make_fullband_downsampled_cube(spw, output_file=None, zarr_path=None,
                                   pixscale=5 * u.arcsec, target_beam=10 * u.arcsec,
                                   parallel=16, block_size=(1, 2048, -1),
                                   zarr_batch_size=64, overwrite=False, verbose=True):
    """Build the full-band, ``pixscale``-pixel, ``target_beam``-beam cube for one spw."""
    filelist = spw_filelist(spw)
    if not filelist:
        raise FileNotFoundError(f"No SPW{spw} feather cubes found")
    weightfilelist = [get_weightfile(fn, spw=spw) for fn in filelist]

    if output_file is None:
        output_file = f'{basepath}/mosaics/cubes/SPW{spw}_fullband_5as_10as_CubeMosaic.fits'
    if zarr_path is None:
        zarr_path = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/SPW{spw}_fullband.zarr'
    if os.path.exists(zarr_path):
        if overwrite:
            shutil.rmtree(zarr_path)
        else:
            raise ValueError(f"zarr_path {zarr_path} exists (pass overwrite=True)")

    header = build_target_header(spw, filelist[0], pixscale=pixscale)
    target_wcs = wcs.WCS(header)
    shape_out = (header['NAXIS3'], header['NAXIS2'], header['NAXIS1'])

    if verbose:
        print(f"SPW{spw}: {len(filelist)} fields -> {shape_out} at {pixscale} pixels", flush=True)

    input_data, input_weights = [], []
    for fn, wfn in zip(filelist, weightfilelist):
        cube = SpectralCube.read(fn, use_dask=True)
        plane = _weight_plane_on_grid(wfn, cube)
        input_data.append((cube.unitless_filled_data[:], cube.wcs))
        input_weights.append(np.broadcast_to(plane[None, :, :], cube.shape))

    array, footprint = reproject_and_coadd(
        input_data, target_wcs, shape_out=shape_out, input_weights=input_weights,
        reproject_function=reproject_interp, combine_function='mean',
        roundtrip_coords=False, parallel=parallel, block_size=block_size,
        return_type='zarr', zarr_path=zarr_path, zarr_batch_size=zarr_batch_size)

    # write the coadd (native-beam, 5") then convolve to the target beam
    if verbose:
        print(f"SPW{spw}: writing coadd -> {output_file}", flush=True)
    fits.PrimaryHDU(data=np.asarray(array), header=header).writeto(output_file, overwrite=overwrite)

    if verbose:
        print(f"SPW{spw}: convolving to {target_beam} beam", flush=True)
    sc = SpectralCube.read(output_file, use_dask=True)
    beam = radio_beam.Beam(target_beam, target_beam, 0 * u.deg)
    sc = sc.convolve_to(beam)
    sc.write(output_file.replace('.fits', f'_{int(target_beam.value)}as.fits'), overwrite=overwrite)
    if verbose:
        print(f"SPW{spw}: done", flush=True)
    return output_file
