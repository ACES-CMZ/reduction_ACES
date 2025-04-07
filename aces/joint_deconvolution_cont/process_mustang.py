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
file = '../data/raw/SgrB2_5pass_1_.0.2_10mJy_10mJy_w_session5_final_smooth4_PlanckCombined_10feb2020.fits'

# opening the FITS file and extracting the first Header Data Unit (HDU)
# HDU is the highest level component of the FITS file structure
hdu = fits.open(file)[0]

# modifying the BMAJ and BMIN header values
# BMAJ and BMIN typically represent the major and minor axes of the beam,
# in degrees. Here they are set to 9 arcseconds.
hdu.header['BMAJ'] = 9/3600
hdu.header['BMIN'] = 9/3600

# setting the BPA (Beam Position Angle) header value to 0
# This angle is usually measured in degrees from north through east
hdu.header['BPA'] = 0

# setting the BUNIT (brightness unit) to 'Jy/beam'
# This indicates that the data is in units of Janskys per beam
hdu.header['BUNIT'] = 'Jy/beam'

# writing the modified HDU back to a FITS file, with the specified name,
# in the specified location, overwriting the file if it already exists.
hdu.writeto('../data/processed/cont_tp.fits', overwrite=True)
