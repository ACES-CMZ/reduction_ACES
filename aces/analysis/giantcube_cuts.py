import numpy as np
import time
from spectral_cube import SpectralCube
from spectral_cube import BooleanArrayMask
import os
from astropy.io import fits
from scipy import ndimage
from radio_beam import Beam
from astropy import units as u
import dask.array as da
from dask_image import ndmorph, ndmeasure

from dask.diagnostics import ProgressBar
from dask.diagnostics import ResourceProfiler
from dask.distributed import progress
import dafits

from tqdm import tqdm

from aces import conf
from aces.imaging.make_mosaic import makepng

basepath = conf.basepath

cubepath = f'{basepath}/mosaics/cubes/'


def daskarr(x):
    try:
        return da.from_array(x)
    except ValueError:
        return x


def copy_with_progress(src, dst, buffer_size=1024**2):
    if os.path.exists(dst):
        print(f"Destination {dst} exists; skipping copy")
        return
    with open(src, 'rb') as fsrc:
        with open(dst, 'wb') as fdst:
            file_size = os.path.getsize(src)
            with tqdm(total=file_size, unit='B', unit_scale=True, desc=f'Copying {src}') as pbar:
                while True:
                    buf = fsrc.read(buffer_size)
                    if not buf:
                        break
                    fdst.write(buf)
                    pbar.update(len(buf))


def get_noise(cube, noise=None, niter=2, threshold1=5, threshold2=3, verbose=True):
    if hasattr(cube, 'rechunk'):
        howargs = {}
    else:
        howargs = {'how': 'ray', 'progressbar': True}
    if noise is None:
        if verbose:
            print(f"Computing mad_std for iteration 0.  time={time.time():0.1f}")
        noise = cube.mad_std(axis=0, **howargs)
    for ii in range(niter):
        if ii == 0:
            mcube = cube.with_mask((cube < threshold1 * noise) & (cube > -threshold1 * noise))
        else:
            mcube = mcube.with_mask((mcube < threshold2 * noise) & (mcube > -threshold2 * noise))
        if verbose:
            print(f"Computing mad_std for iteration {ii}")
            t0 = time.time()
        noise = mcube.mad_std(axis=0, **howargs)
        if verbose:
            print(f'Iteration {ii} noise: {np.nanmean(noise.value)} npix={mcube.mask.include().sum()}.  dt={time.time() - t0}')

    return noise


def get_prunemask_velo(mask, npix=1):

    if npix < 1:
        npix = 1

    roll1 = np.roll(mask, -1 * npix, axis=0)
    roll2 = np.roll(mask, 1 * npix, axis=0)
    mask_ = mask & roll1 & roll2

    return mask_


def get_prunemask_space(mask, npix=10):

    for kk in tqdm(range(mask.shape[0]), desc='Prunemask'):

        mask_ = mask[kk, :, :]
        ll, jj = ndimage.label(mask_)
        counts = ndimage.histogram(ll, 1, jj + 1, jj)
        good_inds = np.arange(1, jj + 1)[counts >= npix]
        mask[kk, :, :] = np.isin(ll, good_inds)

    return mask


def get_prunemask_space_dask(mask, npix=10):

    masks = []
    for kk in tqdm(range(mask.shape[0]), desc='Prunemask'):

        mask_ = mask[kk, :, :]
        ll, jj = ndmeasure.label(mask_)
        jj = jj.compute()
        if jj > 0:
            counts = ndmeasure.histogram(ll, 1, jj + 1, jj)
            # .compute is required here because the dask type is not interpreted as a boolean index array
            good_inds = np.arange(1, jj + 1)[(counts >= npix).compute()]
            masks.append(da.isin(ll, good_inds))
        else:
            masks.append(da.zeros_like(mask_))

    return da.stack(masks)


def get_noedge_mask(cube, iterations=40):
    """
    Erode the edges on a channel-by-channel basis

    Empirically, there are ~35ish pixels along the edge of fields with visibily inflated noise
    """
    good = cube.mask.include()

    if hasattr(cube, 'rechunk'):
        ndi = ndmorph
    else:
        ndi = ndimage

    struct = ndi.generate_binary_structure(3, 1)
    # make cubic binary structure, but ignore the spectral dimension
    struct[0, :, :] = False
    struct[-1, :, :] = False
    good = ndi.binary_erosion(good, structure=struct, iterations=iterations)

    return good


def get_pruned_mask(cube, noise, threshold1=1.0, threshold2=5.0):
    beam_area_pix = get_beam_area_pix(cube)

    if hasattr(cube, 'rechunk'):
        func = get_prunemask_space_dask
        ndi = ndmorph
        niter = 10 # TODO: figure out if this is enough
    else:
        func = get_prunemask_space
        ndi = ndimage
        niter = -1

    signal_mask_l = (cube > noise * threshold1).include()
    signal_mask_l = func(signal_mask_l, npix=beam_area_pix * 3)

    signal_mask_h = (cube > noise * threshold2).include()
    #signal_mask_h = func(signal_mask_h, npix=beam_area_pix)

    signal_mask_both = ndi.binary_dilation(signal_mask_h, iterations=niter, mask=signal_mask_l)

    return signal_mask_both


def get_beam_area_pix(cube):

    beam = Beam(major=cube.header['BMAJ'] * u.deg,
                minor=cube.header['BMIN'] * u.deg,
                pa=cube.header['BPA'] * u.deg)

    pix = np.absolute(cube.header['CDELT1']) * u.deg
    pix_area = pix**2

    beam_area_pix = beam.sr / (pix_area.to(u.sr))

    return beam_area_pix


def do_pvs(cube, molname, mask=None, mompath=f'{basepath}/mosaics/cubes/moments/',
           howargs={}):

    t0 = time.time()
    print(cube, flush=True)
    if hasattr(cube, 'rechunk'):
        print("Dask graph:\n", cube._data.max().__dask_graph__(), flush=True)

        if (cube.shape[1] > 3000) and cube.shape[2] > 5000:
            print("Rechunking")
            cube = cube.rechunk((-1, 500, 500))
            print("Rechunked")

        print(cube, flush=True)
        print("Dask graph:\n", cube._data.max().__dask_graph__(), flush=True)

    print(f"PV peak intensity.  dt={time.time() - t0}", flush=True)
    pv_max = cube.max(axis=1, **howargs)
    pv_max.write(f"{mompath}/{molname}_CubeMosaic_PV_max.fits", overwrite=True)
    makepng(data=pv_max.value, wcs=pv_max.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_max.png",
            stretch='asinh', min_percent=1, max_percent=99.5)

    print(f"PV mean.  dt={time.time() - t0}")
    pv_mean = cube.mean(axis=1, **howargs)
    pv_mean.write(f"{mompath}/{molname}_CubeMosaic_PV_mean.fits", overwrite=True)
    makepng(data=pv_mean.value, wcs=pv_mean.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_mean.png",
            stretch='asinh', min_percent=1, max_percent=99.5)

    print(f"PV max 2.  dt={time.time() - t0}")
    pv_max = cube.max(axis=2, **howargs)
    pv_max.write(f"{mompath}/{molname}_CubeMosaic_PV_b_max.fits", overwrite=True)
    makepng(data=pv_max.value, wcs=pv_max.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_b_max.png",
            stretch='asinh', min_percent=1, max_percent=99.5)

    print(f"PV mean 2.  dt={time.time() - t0}")
    pv_mean = cube.mean(axis=2, **howargs)
    pv_mean.write(f"{mompath}/{molname}_CubeMosaic_PV_b_mean.fits", overwrite=True)
    makepng(data=pv_mean.value, wcs=pv_mean.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_b_mean.png",
            stretch='asinh', min_percent=1, max_percent=99.5)

    if mask is None:
        print("Loading noise map & calculating mask")
        noise = fits.getdata(f"{mompath}/{molname}_CubeMosaic_noisemap.fits")
        mcube = cube.with_mask(cube > noise * cube.unit)
    else:
        print("Using provided mask")
        mcube = cube.with_mask(mask)

    print(f"PV mean with mask.  dt={time.time() - t0}")
    pv_mean_masked = mcube.mean(axis=1, **howargs)
    pv_mean_masked.write(f"{mompath}/{molname}_CubeMosaic_PV_mean_masked.fits", overwrite=True)
    makepng(data=pv_mean_masked.value, wcs=pv_mean_masked.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_mean_masked.png",
            stretch='asinh', min_percent=1, max_percent=99.5)

    print(f"PV mean with mask 2.  dt={time.time() - t0}")
    pv_mean_masked = mcube.mean(axis=2, **howargs)
    pv_mean_masked.write(f"{mompath}/{molname}_CubeMosaic_PV_b_mean_masked.fits", overwrite=True)
    makepng(data=pv_mean_masked.value, wcs=pv_mean_masked.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_b_mean_masked.png",
            stretch='asinh', min_percent=1, max_percent=99.5)

    # 2.5 sigma cut
    noise_2p5 = fits.getdata(f"{mompath}/{molname}_CubeMosaic_dilated_2p5sig_mask.fits")
    mdcube_2p5 = cube.with_mask(cube > noise_2p5 * cube.unit)

    print(f"PV mean.  dt={time.time() - t0}")
    pv_mean_masked = mdcube_2p5.mean(axis=1, **howargs)
    pv_mean_masked.write(f"{mompath}/{molname}_CubeMosaic_PV_mean_masked_2p5.fits", overwrite=True)
    makepng(data=pv_mean_masked.value, wcs=pv_mean_masked.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_pv_mean_masked_2p5.png",
            stretch='asinh', min_percent=1, max_percent=99.5)

    print(f"PV mean 2.  dt={time.time() - t0}")
    pv_mean_masked = mdcube_2p5.mean(axis=2, **howargs)
    pv_mean_masked.write(f"{mompath}/{molname}_CubeMosaic_PV_b_mean_masked_2p5.fits", overwrite=True)
    makepng(data=pv_mean_masked.value, wcs=pv_mean_masked.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_b_mean_masked_2p5.png",
            stretch='asinh', min_percent=1, max_percent=99.5)


def do_all_stats(cube, molname, mompath=f'{basepath}/mosaics/cubes/moments/',
                 howargs={}):
    t0 = time.time()
    print(cube, flush=True)
    if hasattr(cube, 'rechunk'):
        print("Dask graph:\n", cube._data.max().__dask_graph__(), flush=True)

        if False:
            if (cube.shape[1] > 3000) and cube.shape[2] > 5000:
                print("Rechunking")
                cube = cube.rechunk((-1, 500, 500), save_to_tmp_dir=True)
                print("Rechunked")
                # print("Rechunking to zarr")
                # cube = cube.rechunk((-1, 500, 500), save_to_tmp_dir=True)

            print("Rechunked")
            print(cube, flush=True)
        print("Dask graph:\n", cube._data.max().__dask_graph__(), flush=True)

    print(f"Howargs: {howargs}")
    print(f"max.  dt={time.time() - t0}", flush=True)
    mx = cube.max(axis=0, **howargs)
    print(f"Done computing max, now writing. dt={time.time() - t0}")
    mx.write(f"{mompath}/{molname}_CubeMosaic_max.fits", overwrite=True)
    print(f"Done writing max, now pnging. dt={time.time() - t0}")
    makepng(data=mx.value, wcs=mx.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_max.png",
            stretch='asinh', min_percent=0.1, max_percent=99.9)
    header = mx.header
    wcs = mx.wcs
    del mx

    print(f"Noisemap. dt={time.time() - t0}")

    if hasattr(cube, 'rechunk'):
        print(f"madstd (dask).  dt={time.time() - t0}")
        noise_madstd = cube.mad_std(axis=0)
    else:
        print(f"madstd (numpy).  dt={time.time() - t0}")
        noise_madstd = cube.mad_std(axis=0, how='ray', progressbar=True)
    print(f"Done with madstd. Writing noisemap to {mompath}/{molname}_CubeMosaic_madstd.fits.  dt={time.time() - t0}")
    noise_madstd.write(f"{mompath}/{molname}_CubeMosaic_madstd.fits", overwrite=True)
    makepng(data=noise_madstd.value, wcs=noise_madstd.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_madstd.png",
            stretch='asinh', min_percent=0.5, max_percent=99.5)

    noise = get_noise(cube, noise=noise_madstd)
    # mcube = cube.with_mask(cube > noise)

    print(f"Done with noisemap. Writing noisemap to {mompath}/{molname}_CubeMosaic_noisemap.fits.  dt={time.time() - t0}")
    noise.write(f"{mompath}/{molname}_CubeMosaic_noisemap.fits", overwrite=True)
    print(f"Done writing noisemap. dt={time.time() - t0}")
    makepng(data=noise.value, wcs=noise.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_noisemap.png",
            stretch='asinh', min_percent=0.5, max_percent=99.5)

    print(f"masked mom0.  dt={time.time() - t0}")
    mcube = cube.with_mask(cube > noise)
    mom0 = mcube.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)

    print(f"Pruned mask. dt={time.time() - t0}")
    signal_mask_both = get_pruned_mask(cube, noise, threshold1=1.5, threshold2=7.0)
    print(f"Done creating pruned mask. dt={time.time() - t0}")
    if hasattr(signal_mask_both, 'compute'):
        dafits.write(f"{mompath}/{molname}_CubeMosaic_signal_mask_pruned.fits",
                     data=signal_mask_both.astype('int'),
                     header=cube.header,
                     overwrite=True)
    else:
        hdu = fits.PrimaryHDU(data=signal_mask_both.astype('int'), header=header)
        hdu.writeto(f"{mompath}/{molname}_CubeMosaic_signal_mask_pruned.fits", overwrite=True)
        del hdu

    print(f"Dilated mask high-to-low sigma moment 0. dt={time.time() - t0}")
    mdcube_both = cube.with_mask(BooleanArrayMask(signal_mask_both, cube.wcs))
    mom0 = mdcube_both.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_hlsig_dilated_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_hlsig_dilated_masked_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)

    print(f"Dilated mask high-to-low sigma moment 1. dt={time.time() - t0}")
    mom1 = mdcube_both.moment1(axis=0, **howargs)
    mom1.write(f"{mompath}/{molname}_CubeMosaic_masked_hlsig_dilated_mom1.fits", overwrite=True)
    makepng(data=mom1.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_hlsig_dilated_masked_mom1.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)

    print(f"Third dilated mask (5-sigma). dt={time.time() - t0}")
    signal_mask_5p0 = cube > noise * 5.0
    signal_mask_5p0 = ndmorph.binary_dilation(daskarr(signal_mask_5p0.include()), structure=np.ones([1, 3, 3]), iterations=1)

    dilated_mask_space_5p0 = signal_mask_5p0.sum(axis=0)
    fits.PrimaryHDU(data=dilated_mask_space_5p0,
                    header=cube.wcs.celestial.to_header()).writeto(f"{mompath}/{molname}_CubeMosaic_dilated_5p0sig_mask.fits", overwrite=True)

    print(f"Dilated mask 5 sigma moment 0. dt={time.time() - t0}")
    mdcube_5p0 = cube.with_mask(BooleanArrayMask(mask=signal_mask_5p0, wcs=cube.wcs))
    mom0 = mdcube_5p0.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_5p0sig_dilated_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_5p0sig_dilated_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)

    print(f"Computing first signal mask (1-sigma). dt={time.time() - t0}")
    signal_mask = cube > noise

    signal_mask = ndmorph.binary_dilation(daskarr(signal_mask.include()), structure=np.ones([1, 3, 3]), iterations=1)

    dilated_mask_space = signal_mask.sum(axis=0)
    fits.PrimaryHDU(data=dilated_mask_space,
                    header=cube.wcs.celestial.to_header()).writeto(f"{mompath}/{molname}_CubeMosaic_dilated_mask.fits", overwrite=True)

    print(f"Dilated mask moment 0. dt={time.time() - t0}")
    mdcube = cube.with_mask(BooleanArrayMask(mask=signal_mask, wcs=cube.wcs))
    mom0 = mdcube.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_dilated_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_dilated_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)

    print(f"Second dilated mask (2.5-sigma). dt={time.time() - t0}")
    signal_mask_2p5 = cube > noise * 2.5
    signal_mask_2p5 = ndmorph.binary_dilation(daskarr(signal_mask_2p5.include()), structure=np.ones([1, 3, 3]), iterations=1)

    dilated_mask_space_2p5 = signal_mask_2p5.sum(axis=0)
    fits.PrimaryHDU(data=dilated_mask_space_2p5,
                    header=cube.wcs.celestial.to_header()).writeto(f"{mompath}/{molname}_CubeMosaic_dilated_2p5sig_mask.fits", overwrite=True)

    print(f"Dilated mask 2.5 sigma moment 0. dt={time.time() - t0}")
    mdcube_2p5 = cube.with_mask(BooleanArrayMask(mask=signal_mask_2p5, wcs=cube.wcs))
    mom0 = mdcube_2p5.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_2p5sig_dilated_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_2p5sig_dilated_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)

    if hasattr(cube, 'rechunk'):
        print(f"argmax (dask).  dt={time.time() - t0}")
        argmx = cube.argmax(axis=0)
    else:
        print(f"argmax (numpy).  dt={time.time() - t0}")
        argmx = cube.argmax(axis=0, how='ray', progressbar=True)
    print(f"Done computing argmax, now computing vmax. dt={time.time() - t0}")
    vmax = cube.spectral_axis[argmx]
    hdu = fits.PrimaryHDU(data=vmax.value, header=header)
    print(f"Done computing vmax, now writing. dt={time.time() - t0}")
    hdu.writeto(f"{mompath}/{molname}_CubeMosaic_vpeak.fits", overwrite=True)
    print(f"Done writing vpeak, now pnging. dt={time.time() - t0}")
    makepng(data=vmax.value, wcs=wcs, imfn=f"{mompath}/{molname}_CubeMosaic_vpeak.png",
            stretch='asinh', min_percent=0.1, max_percent=99.9)
    del vmax

    print(f"mom0.  dt={time.time() - t0}", flush=True)
    mom0 = cube.moment0(axis=0, **howargs)
    print(f"Done computing mom0, now writing. dt={time.time() - t0}")
    mom0.write(f"{mompath}/{molname}_CubeMosaic_mom0.fits", overwrite=True)
    print(f"Done writing mom0, now pnging. dt={time.time() - t0}")
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)
    del mom0

    print(f"Computing no-edge mask. dt={time.time() - t0}")
    noedge = get_noedge_mask(cube)

    print(f"no-edge max.  dt={time.time() - t0}", flush=True)
    nemx = cube.with_mask(noedge).max(axis=0, **howargs)
    print(f"Done computing no-edge max, now writing. dt={time.time() - t0}")
    nemx.write(f"{mompath}/{molname}_CubeMosaic_edgelessmax.fits", overwrite=True)
    print(f"Done writing no-edge max, now pnging. dt={time.time() - t0}")
    makepng(data=nemx.value, wcs=nemx.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_edgelessmax.png",
            stretch='asinh', min_percent=0.1, max_percent=99.9)

    print(f"No-edge Noisemap. dt={time.time() - t0}")

    if hasattr(cube, 'rechunk'):
        print(f"no-edge madstd (dask).  dt={time.time() - t0}")
        noise_madstd = cube.with_mask(noedge).mad_std(axis=0)
    else:
        print(f"no-edge madstd (numpy).  dt={time.time() - t0}")
        noise_madstd = cube.with_mask(noedge).mad_std(axis=0, how='ray', progressbar=True)
    print(f"Done with no-edge madstd. Writing noisemap to {mompath}/{molname}_CubeMosaic_madstd.fits.  dt={time.time() - t0}")
    noise_madstd.write(f"{mompath}/{molname}_CubeMosaic_edgelessmadstd.fits", overwrite=True)
    makepng(data=noise_madstd.value, wcs=noise_madstd.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_edgelessmadstd.png",
            stretch='asinh', min_percent=0.5, max_percent=99.5)

    return BooleanArrayMask(signal_mask_both, cube.wcs)


def main():
    dodask = os.getenv('USE_DASK')
    if dodask and dodask.lower() == 'false':
        dodask = False
    dods = os.getenv('DOWNSAMPLE')
    if dods and dods.lower() == 'false':
        dods = False
    do_pv = os.getenv('DO_PV')
    if do_pv and do_pv.lower() == 'false':
        do_pv = False

    if os.getenv('MOLNAME'):
        molname = os.getenv('MOLNAME')
    else:
        molname = 'CS21'

    print("Relevant environment variables: "
          f"MOLNAME={os.getenv('MOLNAME')} USE_DASK={os.getenv('USE_DASK')} DOWNSAMPLE={os.getenv('DOWNSAMPLE')} DO_PV={os.getenv('DO_PV')}"
          f"DASK_CLIENT={os.getenv('DASK_CLIENT')} SLURM_MEM_PER_NODE={os.getenv('SLURM_MEM_PER_NODE')} SLURM_NTASKS_PER_NODE={os.getenv('SLURM_NTASKS_PER_NODE')} SLURM_NTASKS={os.getenv('SLURM_NTASKS')} SLURM_STEP_NUM_TASKS={os.getenv('SLURM_STEP_NUM_TASKS')} SLURM_TASKS_PER_NODE={os.getenv('SLURM_TASKS_PER_NODE')}"
    )

    print(f"giantcube_cuts main parameters: molname{molname} dodask={dodask} dods={dods} do_pv={do_pv}")

    cubefilename = f'{cubepath}/{molname}_CubeMosaic_spectrally.fits'
    if not os.path.exists(cubefilename):
        cubefilename = f'{cubepath}/{molname}_CubeMosaic.fits'
    print(f"Using cube {cubefilename}")
    if os.getenv("USE_LOCAL"):
        src_cubefilename = cubefilename
        if os.path.getsize(src_cubefilename) > 3e11:
            # too big for local
            print("Using local (copying to /red)")
            cubefilename = os.path.join('/red/adamginsburg/ACES/workdir/', os.path.basename(src_cubefilename))
        else:
            print("Using local (copying to /tmp)")
            cubefilename = os.path.join(os.getenv('TMPDIR'), os.path.basename(src_cubefilename))
        copy_with_progress(src_cubefilename, cubefilename)

    if not dodask:
        print("Before imports.  Slice-wise reduction", flush=True)

        t0 = time.time()

        cube = SpectralCube.read(cubefilename)
        if cube.shape[0] > 500 and 'H40a' not in cubefilename:
            cube = SpectralCube.read(cubefilename.replace(".fits", "_spectrally.fits"))

        howargs = {'how': 'slice', 'progressbar': True}
        print(f"Non-Dask: how are we computing? {howargs}")
    else:
        howargs = {}
        print("Before imports (using dask)", flush=True)
        for keyword in ("SLURM_NTASKS_PER_NODE", "SLURM_NTASKS", "SLURM_STEP_NUM_TASKS", "SLURM_TASKS_PER_NODE",):
            nworkers = os.getenv(keyword)
            if nworkers is not None:
                nworkers = int(nworkers)
                break
        if nworkers is None:
            nworkers = 1
        else:
            nworkers = int(nworkers)
        print(f"Using {nworkers} workers", flush=True)

        memory_limit = os.getenv('SLURM_MEM_PER_NODE')
        if memory_limit:
            memory_limit = f"{int(int(memory_limit) / nworkers)}MB"
        print(f"Memory limit: {memory_limit}")

        cube = SpectralCube.read(cubefilename, use_dask=True)
        if cube.shape[0] > 500 and 'H40a' not in cubefilename:
            cube = SpectralCube.read(cubefilename.replace(".fits", "_spectrally.fits"), use_dask=True)

        t0 = time.time()

        import dask
        threads = os.getenv('DASK_CLIENT') == 'threads'
        # dask_does_slurm supercedes threads
        dask_does_slurm = os.getenv('DASK_CLIENT') == 'slurm'

        from dask.config import set as dask_set

        dask_set({"array.chunk-size": "1GB"})

        if dask_does_slurm:
            from dask_jobqueue import SLURMCluster
            from dask.distributed import Client
            cluster = SLURMCluster(
                job_extra_directives=['--qos="astronomy-dept-b"'],
                account="astronomy-dept",
                cores=16,
                memory="64 GB"
            )
            cluster.scale(jobs=10)  # ask for 10 jobs
            client = Client(cluster)
            scheduler = cube.use_dask_scheduler(client, num_workers=nworkers)
            print("Dask-Slurm")
            print(f"dashboard: {cluster.dashboard_link}", flush=True)
            print(f"Using dask scheduler {client} ({'threads' if threads else 'processes'})", flush=True)
            print(f"Client info: {client.scheduler_info()}", flush=True)
            print(f"Dask client number of workers: {len(client.scheduler_info()['workers'])}")
            print(f"Using scheduler {scheduler}", flush=True)
        elif threads:
            print("Using threaded scheduler (no dask.distributed LocalCluster)")
            scheduler = cube.use_dask_scheduler('threads', num_workers=nworkers)
            dask.config.set(scheduler='threads')
            pbar = ProgressBar()
            pbar.register()
        else:
            from dask.distributed import Client, LocalCluster
            threads_per_worker = nworkers if threads else 1
            cluster = LocalCluster(n_workers=nworkers, memory_limit=memory_limit,
                                   processes=not threads, threads_per_worker=threads_per_worker)
            client = Client(cluster)
            print(f"dashboard: {cluster.dashboard_link}", flush=True)

            # switch to our LocalCluster scheduler
            scheduler = cube.use_dask_scheduler(client, num_workers=nworkers)
            print(f"Using dask scheduler {client} ({'threads' if threads else 'processes'})", flush=True)
            print(f"Client info: {client.scheduler_info()}", flush=True)
            print(f"Dask client number of workers: {len(client.scheduler_info()['workers'])}")
            print(f"Using scheduler {scheduler}", flush=True)

    signal_mask_both = do_all_stats(cube, molname=molname, howargs=howargs)

    if do_pv:
        do_pvs(cube, molname=molname, howargs=howargs, mask=signal_mask_both)

    if dods:
        print(f"Downsampling.  dt={time.time() - t0}")
        from aces.imaging.make_mosaic import make_downsampled_cube, basepath
        make_downsampled_cube(f'{cubepath}/{molname}_CubeMosaic.fits', f'{cubepath}/{molname}_CubeMosaic_downsampled9.fits',
                              smooth_beam=12 * u.arcsec)


if __name__ == "__main__":
    for key in os.environ:
        if 'SLURM' in key:
            print(f"{key} = {os.getenv(key)}")
    print("\n\n")
    main()
