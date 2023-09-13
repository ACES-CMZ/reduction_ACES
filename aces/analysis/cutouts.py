from spectral_cube import SpectralCube
import numpy as np
from astropy.visualization import quantity_support
from astropy import units as u
from astropy import wcs
import matplotlib.pyplot as plt
from astropy.utils import data
from reproject import reproject_exact
from astropy.io import fits
import regions
import glob, os

reg = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions/nircam_cloudc_fov.reg')

for fn in glob.glob("/orange/adamginsburg/ACES/mosaics/cubes/*CubeMosaic.fits"):
    cube = SpectralCube.read(fn)
    scube = cube.subcube_from_regions(reg, minimize=False)
    print(fn)
    scube.write(f'/orange/adamginsburg/ACES/20230918Meeting/cutouts/cloudc_{os.path.basename(fn).replace("_CubeMosaic","_cutout")}')



