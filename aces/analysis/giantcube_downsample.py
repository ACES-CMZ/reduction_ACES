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
from aces.imaging.make_mosaic import make_downsampled_cube, basepath

basepath = conf.basepath

cubepath = f'{basepath}/mosaics/cubes/'

if os.getenv('MOLNAME'):
    molname = os.getenv('MOLNAME')
else:
    molname = 'CS21'


print(f"Downsampling {molname}.")
make_downsampled_cube(f'{cubepath}/{molname}_CubeMosaic.fits',
                      f'{cubepath}/{molname}_CubeMosaic_downsampled9.fits',
                      smooth_beam=12*u.arcsec,
                      overwrite=True,
                      use_dask=False
                      )
