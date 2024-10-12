#from dask.diagnostics import ProgressBar
#pbar = ProgressBar(minimum=20) # don't show pbar <20s
#pbar.register()
import numpy as np
import time
from spectral_cube import SpectralCube
from spectral_cube import BooleanArrayMask
import os
from astropy.io import fits
from scipy import ndimage
from radio_beam import Beam
from astropy import units as u

import shutil
import os
from tqdm import tqdm

from aces import conf
from aces.imaging.make_mosaic import makepng

basepath = conf.basepath

cubepath = f'{basepath}/mosaics/cubes/'


def copy_with_progress(src, dst, buffer_size=1024*1024):
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


def get_noise(cube, niter=5, threshold1=3, threshold2=3, verbose=True):
    if verbose:
        print(f"Computing mad_std for iteration 0.  time={time.time():0.1f}")
    noise = cube.mad_std(axis=0)
    for ii in range(niter):
        if ii == 0:
            mcube = cube.with_mask(cube < threshold1 * noise)
        else:
            mcube = mcube.with_mask(mcube < threshold2 * noise)
        if verbose:
            print(f"Computing mad_std for iteration {ii}")
            t0 = time.time()
        noise = mcube.mad_std(axis=0)
        if verbose:
            print(f'Iteration {ii} noise: {np.nanmean(noise.value)}.  dt={time.time() - t0}')

    return noise


def get_prunemask_velo(mask, npix=1):

    if npix < 1:
        npix = 1

    roll1 = np.roll(mask, -1 * npix, axis=0)
    roll2 = np.roll(mask, 1 * npix, axis=0)
    mask_ = mask & roll1 & roll2

    return mask_


def get_prunemask_space(mask, npix=10):

    for kk in tqdm.tqdm(range(mask.shape[0], desc='Prunemask')):

        mask_ = mask[kk, :, :]
        ll, jj = ndimage.label(mask_)
        hist = ndimage.measurements.histogram(ll, 0, jj+1, jj+1)
        os = ndimage.find_objects(ll)

        for ii in range(jj):
            if hist[ii + 1] < npix:
                mask_[os[ii]] = 0

        mask[kk, :, :] = mask_

    return mask


def get_beam_area_pix(cube):

    beam = Beam(major=cube.header['BMAJ'] * u.deg,
                    minor=cube.header['BMIN'] * u.deg,
                    pa=cube.header['BPA'] * u.deg)

    pix = np.absolute(cube.header['CDELT1']) * u.deg
    pix_area = pix**2

    beam_area_pix = beam.sr / (pix_area.to(u.sr))

    return beam_area_pix


def do_all_stats(cube, molname, mompath=f'{basepath}/mosaics/cubes/moments/',
                 dopv=True, dods=True, howargs={}):
    t0 = time.time()
    print(cube, flush=True)
    if hasattr(cube, 'rechunk'):
        print("Dask graph:\n", cube._data.max().__dask_graph__(), flush=True)

        cube = cube.rechunk((-1, 1000, 1000), save_to_tmp_dir=True)

        print("Rechunked")
        print(cube, flush=True)
        print("Dask graph:\n", cube._data.max().__dask_graph__(), flush=True)

    print(f"mom0.  dt={time.time() - t0}", flush=True)
    mom0 = cube.moment0(axis=0, **howargs)
    print(f"Done computing mom0, now writing. dt={time.time() - t0}")
    mom0.write(f"{mompath}/{molname}_CubeMosaic_mom0.fits", overwrite=True)
    print(f"Done writing mom0, now pnging. dt={time.time() - t0}")
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)
    del mom0


    print(f"max.  dt={time.time() - t0}", flush=True)
    mx = cube.max(axis=0, **howargs)
    print(f"Done computing max, now writing. dt={time.time() - t0}")
    mx.write(f"{mompath}/{molname}_CubeMosaic_max.fits", overwrite=True)
    print(f"Done writing max, now pnging. dt={time.time() - t0}")
    makepng(data=mx.value, wcs=mx.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_max.png",
            stretch='asinh', min_percent=0.1, max_percent=99.9)
    del mx

    print(f"argmax.  dt={time.time() - t0}")
    if hasattr(cube, 'rechunk'):
        argmx = cube.argmax(axis=0)
    else:
        argmx = cube.argmax(axis=0, how='ray', progressbar=True)
    print(f"Done computing argmax, now computing vmax. dt={time.time() - t0}")
    vmax = cube.spectral_axis[argmx]
    hdu = mx.hdu
    hdu.data = vmax.value
    print(f"Done computing vmax, now writing. dt={time.time() - t0}")
    hdu.writeto(f"{mompath}/{molname}_CubeMosaic_vpeak.fits", overwrite=True)
    # use mx.wcs
    print(f"Done writing vpeak, now pnging. dt={time.time() - t0}")
    makepng(data=vmax.value, wcs=mx.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_vpeak.png",
            stretch='asinh', min_percent=0.1, max_percent=99.9)

    if dopv:
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

    if dods:
        print(f"Downsampling.  dt={time.time() - t0}")
        from aces.imaging.make_mosaic import make_downsampled_cube, basepath
        make_downsampled_cube(f'{cubepath}/{molname}_CubeMosaic.fits', f'{cubepath}/{molname}_CubeMosaic_downsampled9.fits',
                              smooth_beam=12*u.arcsec)

    print(f"Noisemap. dt={time.time() - t0}")

    noise = get_noise(cube)
    # noise = cube.mad_std(axis=0)
    # mcube = cube.with_mask(cube > noise)

    print(f"Done with noisemap. Writing noisemap to {mompath}/{molname}_CubeMosaic_madstd.fits.  dt={time.time() - t0}")
    noise.write(f"{mompath}/{molname}_CubeMosaic_madstd.fits", overwrite=True)
    print("Done writing noisemap. dt={time.time() - t0}")
    makepng(data=noise.value, wcs=noise.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_madstd.png",
            stretch='asinh', min_percent=0.5, max_percent=99.5)

    print(f"masked mom0.  dt={time.time() - t0}")
    #try:
    #    std = cube.mad_std()
    #except ValueError:
    #    # mad_std requires whole cube in memory; we can't afford that
    #    # instead, do a cheap version of sigma clipping
    #    std = cube.std()
    #    std = cube.with_mask(cube < std * 5).std()
    mcube = cube.with_mask(cube > noise)
    mom0 = mcube.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)

    if dopv:
        print(f"PV mean.  dt={time.time() - t0}")
        pv_mean_masked = mcube.mean(axis=1, **howargs)
        pv_mean_masked.write(f"{mompath}/{molname}_CubeMosaic_PV_mean_masked.fits", overwrite=True)
        makepng(data=pv_mean.value, wcs=pv_mean.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_mean_masked.png",
                stretch='asinh', min_percent=1, max_percent=99.5)

        print(f"PV mean 2.  dt={time.time() - t0}")
        pv_mean_masked = mcube.mean(axis=2, **howargs)
        pv_mean_masked.write(f"{mompath}/{molname}_CubeMosaic_PV_b_mean_masked.fits", overwrite=True)
        makepng(data=pv_mean.value, wcs=pv_mean.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_b_mean_masked.png",
                stretch='asinh', min_percent=1, max_percent=99.5)

    # do last b/c it doens't work right now
    """
    Traceback (most recent call last):
      File "/orange/adamginsburg/ACES/reduction_ACES/aces/analysis/giantcube_cuts.py", line 123, in <module>
        signal_mask = ndmorph.binary_dilation(signal_mask, structure=np.ones([3, 3, 3]), iterations=1)
      File "/orange/adamginsburg/miniconda3/envs/python39/lib/python3.9/site-packages/dask_image/ndmorph/__init__.py", line 57, in binary_dilation
        dispatch_binary_dilation(image),
      File "/orange/adamginsburg/miniconda3/envs/python39/lib/python3.9/site-packages/dask_image/dispatch/_dispatcher.py", line 23, in __call__
        meth = self.dispatch(datatype)
      File "/orange/adamginsburg/miniconda3/envs/python39/lib/python3.9/site-packages/dask/utils.py", line 635, in dispatch
        raise TypeError(f"No dispatch for {cls}")
    TypeError: No dispatch for <class 'spectral_cube.masks.LazyComparisonMask'>
    """
    from dask_image import ndmorph
    print(f"Computing first signal mask (1-sigma). dt={time.time() - t0}")
    signal_mask = cube > noise

    print(f"Third dilated mask (5-sigma). dt={time.time() - t0}")
    signal_mask_5p0 = cube > noise * 5.0
    signal_mask_5p0 = ndmorph.binary_dilation(signal_mask_5p0.include(), structure=np.ones([1, 3, 3]), iterations=1)

    dilated_mask_space_5p0 = signal_mask_5p0.sum(axis=0)
    fits.PrimaryHDU(data=dilated_mask_space_5p0,
                    header=cube.wcs.celestial.to_header()).writeto(f"{mompath}/{molname}_CubeMosaic_dilated_5p0sig_mask.fits", overwrite=True)

    print(f"Dilated mask 5 sigma moment 0. dt={time.time() - t0}")
    mdcube_5p0 = cube.with_mask(BooleanArrayMask(mask=signal_mask_5p0, wcs=cube.wcs))
    mom0 = mdcube_5p0.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_5p0sig_dilated_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_5p0sig_dilated_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)



    signal_mask = ndmorph.binary_dilation(signal_mask.include(), structure=np.ones([1, 3, 3]), iterations=1)

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
    signal_mask_2p5 = ndmorph.binary_dilation(signal_mask_2p5.include(), structure=np.ones([1, 3, 3]), iterations=1)

    dilated_mask_space_2p5 = signal_mask_2p5.sum(axis=0)
    fits.PrimaryHDU(data=dilated_mask_space_2p5,
                    header=cube.wcs.celestial.to_header()).writeto(f"{mompath}/{molname}_CubeMosaic_dilated_2p5sig_mask.fits", overwrite=True)


    print(f"Dilated mask 2.5 sigma moment 0. dt={time.time() - t0}")
    mdcube_2p5 = cube.with_mask(BooleanArrayMask(mask=signal_mask_2p5, wcs=cube.wcs))
    mom0 = mdcube_2p5.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_2p5sig_dilated_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_2p5sig_dilated_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)


    print(f"Third dilated mask (5-sigma). dt={time.time() - t0}")
    signal_mask_5p0 = cube > noise * 5.0
    signal_mask_5p0 = ndmorph.binary_dilation(signal_mask_5p0.include(), structure=np.ones([1, 3, 3]), iterations=1)

    dilated_mask_space_5p0 = signal_mask_5p0.sum(axis=0)
    fits.PrimaryHDU(data=dilated_mask_space_5p0,
                    header=cube.wcs.celestial.to_header()).writeto(f"{mompath}/{molname}_CubeMosaic_dilated_5p0sig_mask.fits", overwrite=True)

    print(f"Dilated mask 5 sigma moment 0. dt={time.time() - t0}")
    mdcube_5p0 = cube.with_mask(BooleanArrayMask(mask=signal_mask_5p0, wcs=cube.wcs))
    mom0 = mdcube_5p0.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_5p0sig_dilated_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_5p0sig_dilated_mom0.png",
            stretch='asinh', vmin=-0.1, max_percent=99.5)

    print(f"Fourth mask (high-to-low sigma). dt={time.time() - t0}")
    signal_mask_l = cube > noise * 1.5
    signal_mask_h = cube > noise * 3.0

    # Not super obvious we want to include this for all cases...
    # As for the 3km/s having a +-1 connected pixels in velocity may be too aggressive..
    # signal_mask_l = get_prunemask_velo(signal_mask_l.include(), npix=1)
    # signal_mask_h = get_prunemask_velo(signal_mask_h.include(), npix=1)

    # signal_mask_l = get_prunemask_space(signal_mask_l.include(), npix=5)
    # signal_mask_h = get_prunemask_space(signal_mask_h.include(), npix=5)
    print(f"Pruned mask. dt={time.time() - t0}")
    beam_area_pix = get_beam_area_pix(cube) # Remove small islands of noise
    signal_mask_l = get_prunemask_space(signal_mask_l.include(), npix=beam_area_pix * 3)
    signal_mask_h = get_prunemask_space(signal_mask_h.include(), npix=beam_area_pix * 3)

    signal_mask_both = ndmorph.binary_dilation(signal_mask_h, iterations=-1, mask=signal_mask_l)

    print(f"Dilated mask high-to-low sigma moment 0. dt={time.time() - t0}")
    mdcube_both = cube.with_mask(signal_mask_both)
    mom0 = mdcube_both.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_hlsig_dilated_mom0.fits", overwrite=True)

    if dopv:
        print(f"PV mean.  dt={time.time() - t0}")
        pv_mean_masked = mdcube_2p5.mean(axis=1, **howargs)
        pv_mean_masked.write(f"{mompath}/{molname}_CubeMosaic_PV_mean_masked_2p5.fits", overwrite=True)
        makepng(data=pv_mean.value, wcs=pv_mean.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_mean_masked_2p5.png",
                stretch='asinh', min_percent=1, max_percent=99.5)

        print(f"PV mean 2.  dt={time.time() - t0}")
        pv_mean_masked = mdcube_2p5.mean(axis=2, **howargs)
        pv_mean_masked.write(f"{mompath}/{molname}_CubeMosaic_PV_b_mean_masked_2p5.fits", overwrite=True)
        makepng(data=pv_mean.value, wcs=pv_mean.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_PV_b_mean_masked_2p5.png",
                stretch='asinh', min_percent=1, max_percent=99.5)


def main():
    dodask = os.getenv('USE_DASK')
    if dodask and dodask.lower() == 'false':
        dodask = False
    dods = os.getenv('DOWNSAMPLE')
    dopv = os.getenv('DO_PV')

    if os.getenv('MOLNAME'):
        molname = os.getenv('MOLNAME')
    else:
        molname = 'CS21'

    cubefilename = f'{cubepath}/{molname}_CubeMosaic.fits'
    if os.getenv("USE_LOCAL"):
        src_cubefilename = cubefilename
        cubefilename = os.path.join(os.getenv('TMPDIR'), os.path.basename(cubefilename))
        copy_with_progress(src_cubefilename, cubefilename)

    if not dodask:
        print("Before imports.  Slice-wise reduction", flush=True)

        t0 = time.time()

        cube = SpectralCube.read(cubefilename)

        howargs = {'how': 'slice', 'progressbar': True}
        print(f"Non-Dask: how are we computing? {howargs}")
    else:
        print("Before imports (using dask)", flush=True)
        nworkers = os.getenv('SLURM_NTASKS_PER_NODE')
        if nworkers is None:
            nworkers = 1
        else:
            nworkers = int(nworkers)
        print(f"Using {nworkers} workers", flush=True)

        memory_limit = os.getenv('SLURM_MEM_PER_CPU')
        if memory_limit:
            memory_limit = f"{memory_limit}MB"
        print(f"Memory limit: {memory_limit}")

        cube = SpectralCube.read(cubefilename, use_dask=True)

        t0 = time.time()

        import dask
        #dask.config.set(scheduler='threads')
        from dask.distributed import Client, LocalCluster
        threads = os.getenv('DASK_CLIENT') == 'threads'
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

        howargs = {}

        from dask.diagnostics import ProgressBar
        from dask.diagnostics import ResourceProfiler
        from dask.distributed import progress

    do_all_stats(cube, molname=molname, dopv=dopv, dods=dods, howargs=howargs)


if __name__ == "__main__":
    main()