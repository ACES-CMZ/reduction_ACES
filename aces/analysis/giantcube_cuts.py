#from dask.diagnostics import ProgressBar
#pbar = ProgressBar(minimum=20) # don't show pbar <20s
#pbar.register()
import time
from spectral_cube import SpectralCube
import os

from aces import conf

basepath = conf.basepath

cubepath = f'{basepath}/mosaics/cubes/'
mompath = f'{basepath}/mosaics/cubes/moments/'

if __name__ == "__main__":
    dodask = os.getenv('USE_DASK')

    if os.getenv('MOLNAME'):
        molname = os.getenv('MOLNAME')
    else:
        molname = 'CS21'

    if not dodask:
        print("Before imports.  Slice-wise reduction", flush=True)

        t0 = time.time()

        cube = SpectralCube.read(f'{cubepath}/{molname}_CubeMosaic.fits')
        print(cube)

        print(f"PV peak intensity.  dt={time.time()-t0}", flush=True)
        pv_max = cube.max(axis=1, how='slice')
        pv_max.write(f"{mompath}/{molname}_CubeMosaic_PV_max.fits", overwrite=True)

        print(f"PV mean.  dt={time.time()-t0}")
        pv_mean = cube.mean(axis=1, how='slice')
        pv_mean.write(f"{mompath}/{molname}_CubeMosaic_PV_mean.fits", overwrite=True)

        print(f"mom0.  dt={time.time()-t0}")
        mom0 = cube.moment0(axis=0, how='slice')
        mom0.write(f"{mompath}/{molname}_CubeMosaic_mom0.fits", overwrite=True)

        print(f"max.  dt={time.time()-t0}")
        mx = cube.max(axis=0, how='slice')
        mx.write(f"{mompath}/{molname}_CubeMosaic_max.fits", overwrite=True)

        print(f"PV max 2.  dt={time.time()-t0}")
        pv_max = cube.max(axis=2, how='slice')
        pv_max.write(f"{mompath}/{molname}_CubeMosaic_PV_b_max.fits", overwrite=True)

        print(f"PV mean 2.  dt={time.time()-t0}")
        pv_mean = cube.mean(axis=2, how='slice')
        pv_mean.write(f"{mompath}/{molname}_CubeMosaic_PV_b_mean.fits", overwrite=True)

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
        print(cube)

        print(f"PV peak intensity.  dt={time.time()-t0}", flush=True)
        pv_max = cube.max(axis=1)
        pv_max.write(f"{mompath}/{molname}_CubeMosaic_PV_max.fits", overwrite=True)

        print(f"PV mean.  dt={time.time()-t0}")
        pv_mean = cube.mean(axis=1)
        pv_mean.write(f"{mompath}/{molname}_CubeMosaic_PV_mean.fits", overwrite=True)

        print(f"mom0.  dt={time.time()-t0}")
        mom0 = cube.moment0(axis=0)
        mom0.write(f"{mompath}/{molname}_CubeMosaic_mom0.fits", overwrite=True)

        print(f"max.  dt={time.time()-t0}")
        mx = cube.max(axis=0)
        mx.write(f"{mompath}/{molname}_CubeMosaic_max.fits", overwrite=True)

        print(f"PV max 2.  dt={time.time()-t0}")
        pv_max = cube.max(axis=2)
        pv_max.write(f"{mompath}/{molname}_CubeMosaic_PV_b_max.fits", overwrite=True)

        print(f"PV mean 2.  dt={time.time()-t0}")
        pv_mean = cube.mean(axis=2)
        pv_mean.write(f"{mompath}/{molname}_CubeMosaic_PV_b_mean.fits", overwrite=True)

    print("Downsampling")
    from aces.imaging.make_mosaic import make_downsampled_cube, basepath
    make_downsampled_cube(f'{cubepath}/{molname}_CubeMosaic.fits', f'{cubepath}/{molname}_CubeMosaic_downsampled9.fits')
