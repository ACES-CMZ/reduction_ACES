# importing required libraries
from glob import glob
from astropy.io import fits
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
import numpy as np
from tqdm.auto import tqdm
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.units as u
import gc
import os

# defining the path of the FITS file to be opened and processed
file = '../data/raw/12m_continuum_commonbeam_circular_reimaged_mosaic.fits'

# opening the FITS file and extracting the first Header Data Unit (HDU)
# HDU is the highest level component of the FITS file structure
hdu = fits.open(file)[0]

# setting the BUNIT (brightness unit) to 'Jy/beam'
# This indicates that the data is in units of Janskys per beam
hdu.header['BUNIT'] = 'Jy/beam'

# writing the modified HDU back to a FITS file, with the specified name,
# in the specified location, overwriting the file if it already exists.
hdu.writeto('../data/processed/cont_12m.fits', overwrite=True)
