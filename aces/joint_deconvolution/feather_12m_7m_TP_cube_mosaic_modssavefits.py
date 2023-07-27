# Import the necessary libraries
import os
import glob
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
from casatasks import imhead, exportfits, imtrans, feather, imreframe, importfits, imregrid, rmtables, imrebin
from spectral_cube import SpectralCube
from tqdm import tqdm 
import astropy.units as u
import gc 
import warnings
warnings.filterwarnings('ignore')


# Function to convert the data in an HDU (Header Data Unit) to float32 format
def convert_to_float32(hdu):
    hdu.data = hdu.data.astype('float32')
    return(hdu)


# Function to crop a FITS cube to a specific velocity range
def cubeconvert_K_kms(ACES_WORKDIR, MOLECULE):

    fits_file =  f'{ACES_WORKDIR}/{MOLECULE}.TP_7M_12M_weighted_mosaic.fits'

    # Load the FITS cube file
    cube = SpectralCube.read(fits_file)
    cube.allow_huge_operations=True

    # Convert the cube to velocity space using the given rest frequency
    cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')
    
    # Crop the cube further to the minimal enclosing subcube
    cube = cube.minimal_subcube()

    # Convert the cube to Kelvin units
    cube = cube.to(u.K)

    # Convert the cube to an HDU
    hdu = cube.hdu

    # Convert the HDU data to float32 format
    hdu = fits.PrimaryHDU(hdu.data, hdu.header)
    hdu = convert_to_float32(hdu) 

    outputfile = fits_file.replace('.fits', '.K.kms.fits')
    print(f"[INFO] Created and saved weighted mosaic to {outputfile}")
    hdu.writeto(outputfile, overwrite=True)

    return None


def rebin(ACES_WORKDIR, MOLECULE, factor=3, overwrite=True):
    """

    """

    input_fits =  f'{ACES_WORKDIR}/{MOLECULE}.TP_7M_12M_weighted_mosaic.K.kms.fits'

    # Define the names of the intermediate images
    input_image = input_fits.replace('.fits', '.tmp.image')
    regrid_image = input_fits.replace('.fits', '.tmp.regrid.image')
    output_fits = input_fits.replace('.fits', '.rebin.fits')

    # Remove any pre-existing intermediate images
    rmtables([input_image, regrid_image])

    # Import the .fits files into CASA images
    importfits(fitsimage=input_fits, imagename=input_image, overwrite=overwrite)

    # Rebin the image by some factor 
    imrebin(imagename=input_image, outfile=regrid_image, factor=[factor,factor,1], overwrite=True)

    # Export the regridded image to a .fits file
    exportfits(imagename=regrid_image, fitsimage=output_fits, overwrite=overwrite, velocity=True)

    # Clean up by removing the intermediate images
    rmtables([input_image, regrid_image])
