#from dask.diagnostics import ProgressBar
#pbar = ProgressBar(minimum=20) # don't show pbar <20s
#pbar.register()
import numpy as np
import time
from spectral_cube import SpectralCube
from spectral_cube import BooleanArrayMask
import os

from aces import conf
from aces.imaging.make_mosaic import makepng

basepath = conf.basepath

cubepath = f'{basepath}/mosaics/cubes/'
mompath = f'{basepath}/mosaics/cubes/moments/'

if __name__ == "__main__":
    dodask = os.getenv('USE_DASK')
    if dodask and dodask.lower() == 'false':
        dodask = False
    dods = os.getenv('DOWNSAMPLE')
    dopv = os.getenv('DO_PV')

    if os.getenv('MOLNAME'):
        molname = os.getenv('MOLNAME')
    else:
        molname = 'CS21'

    if not dodask:
        print("Before imports.  Slice-wise reduction", flush=True)

        t0 = time.time()

        cube = SpectralCube.read(f'{cubepath}/{molname}_CubeMosaic.fits')

        howargs = {'how': 'slice'}
    else:
        print("Before imports (using dask)", flush=True)

        import dask
        #dask.config.set(scheduler='threads')
        from dask.distributed import Client, LocalCluster
        cluster = LocalCluster()
        client = Client(cluster)
        print(f"dashboard: {cluster.dashboard_link}", flush=True)

        t0 = time.time()

        cube = SpectralCube.read(f'{cubepath}/{molname}_CubeMosaic.fits', use_dask=True)

        howargs = {}

    print(cube, flush=True)

    print(f"mom0.  dt={time.time() - t0}", flush=True)
    mom0 = cube.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_mom0.png",
            stretch='asinh', min_cut=-0.1, max_percent=99.5)

    print(f"max.  dt={time.time() - t0}", flush=True)
    mx = cube.max(axis=0, **howargs)
    mx.write(f"{mompath}/{molname}_CubeMosaic_max.fits", overwrite=True)
    makepng(data=mx.value, wcs=mx.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_max.png",
            stretch='asinh', min_percent=0.1, max_percent=99.9)

    print(f"argmax.  dt={time.time() - t0}")
    argmx = cube.argmax(axis=0, **howargs)
    vmax = cube.spectral_axis[argmx]
    hdu = mx.hdu
    hdu.data = vmax.value
    hdu.writeto(f"{mompath}/{molname}_CubeMosaic_vpeak.fits", overwrite=True)
    # use mx.wcs
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
        print("Downsampling")
        from aces.imaging.make_mosaic import make_downsampled_cube, basepath
        make_downsampled_cube(f'{cubepath}/{molname}_CubeMosaic.fits', f'{cubepath}/{molname}_CubeMosaic_downsampled9.fits')

    print(f"Noisemap. dt={time.time() - t0}")
    noise = cube.mad_std(axis=0)
    noise.write(f"{mompath}/{molname}_CubeMosaic_madstd.fits", overwrite=True)
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
            stretch='asinh', min_cut=-0.1, max_percent=99.5)

    mx = mcube.max(axis=0, **howargs)
    mx.write(f"{mompath}/{molname}_CubeMosaic_masked_max.fits", overwrite=True)
    makepng(data=mx.value, wcs=mx.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_max.png",
            stretch='asinh', min_cut=-0.1, max_percent=99.5)

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
    signal_mask = cube > noise
    signal_mask = ndmorph.binary_dilation(signal_mask.include(), structure=np.ones([1, 3, 3]), iterations=1)
    mdcube = cube.with_mask(BooleanArrayMask(mask=signal_mask, wcs=cube.wcs))
    mom0 = mdcube.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_dilated_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_dilated_mom0.png",
            stretch='asinh', min_cut=-0.1, max_percent=99.5)

    mx = mdcube.max(axis=0, **howargs)
    mx.write(f"{mompath}/{molname}_CubeMosaic_masked_dilated_max.fits", overwrite=True)
    makepng(data=mx.value, wcs=mx.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_dilated_max.png",
            stretch='asinh', min_cut=-0.1, max_percent=99.5)

    signal_mask_2p5 = cube > noise * 2.5
    signal_mask_2p5 = ndmorph.binary_dilation(signal_mask_2p5.include(), structure=np.ones([1, 3, 3]), iterations=1)
    mdcube_2p5 = cube.with_mask(BooleanArrayMask(mask=signal_mask_2p5, wcs=cube.wcs))
    mom0 = mdcube_2p5.moment0(axis=0, **howargs)
    mom0.write(f"{mompath}/{molname}_CubeMosaic_masked_2p5sig_dilated_mom0.fits", overwrite=True)
    makepng(data=mom0.value, wcs=mom0.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_2p5sig_dilated_mom0.png",
            stretch='asinh', min_cut=-0.1, max_percent=99.5)

    mx = mdcube_2p5.max(axis=0, **howargs)
    mx.write(f"{mompath}/{molname}_CubeMosaic_masked_2p5sig_dilated_max.fits", overwrite=True)
    makepng(data=mx.value, wcs=mx.wcs, imfn=f"{mompath}/{molname}_CubeMosaic_masked_2p5sig_dilated_max.png",
            stretch='asinh', min_cut=-0.1, max_percent=99.5)

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
