# importing required libraries
from glob import glob
from astropy.io import fits
from reproject import reproject_interp
import numpy as np
from tqdm.auto import tqdm
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
import astropy.units as u
import gc
import os

# defining the path of the FITS file to be opened and processed
file = './../data/feathered/cont_12mtp.fits'

# opening the FITS file and extracting the first Header Data Unit (HDU)
hdu_12mtp = fits.open(file)[0]

# removing the singleton dimensions of the data
hdu_12mtp.data = np.squeeze(hdu_12mtp.data)

# removing some specific header entries
del hdu_12mtp.header['*3*']
del hdu_12mtp.header['*4*']

# writing the modified HDU back to a FITS file, with the specified name
hdu_12mtp.writeto('./../data/feathered/cont_12mtp_final.fits')

# defining the path of another FITS file to be opened
file = './../data/processed/cont_tp_cropped.fits'

# opening the second FITS file and extracting the first HDU
hdu_tp = fits.open(file)[0]

# reprojecting the second HDU to match the coordinate system of the first HDU
# and ignoring the footprint (the second output of the reproject_interp function)
data, _ = reproject_interp(hdu_tp, hdu_12mtp.header)

# calculating the ratio of the beam sizes (based on the major and minor axes)
# of the two HDUs
beam_ratio = (hdu_12mtp.header['BMAJ']*hdu_12mtp.header['BMIN'])/(hdu_tp.header['BMAJ']*hdu_tp.header['BMIN'])

# multiplying the reprojected data by the beam ratio
data = data * beam_ratio

# replacing NaN values in the first HDU's data with corresponding values from the reprojected data
hdu_12mtp.data[np.isnan(hdu_12mtp.data)] = data[np.isnan(hdu_12mtp.data)]

# writing the final modified HDU back to a new FITS file
hdu_12mtp.writeto('./../data/feathered/cont_12mtp_final_filled.fits')
