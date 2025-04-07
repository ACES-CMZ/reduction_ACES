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
from casatasks import imhead, exportfits, imtrans, feather, imreframe, importfits, imregrid, rmtables, imsmooth, rmtables, imrebin
from tqdm.auto import tqdm
import astropy.units as u
import gc 
import warnings
from astropy.convolution import Gaussian2DKernel
from radio_beam import Beam
from astropy import units as u
warnings.filterwarnings('ignore')


def create_fake_hdus(hdus, j):
    '''
    Function to create fake HDUs (Header Data Units) for use in reprojecting.
    
    Inputs:
    - hdus: a list of HDU objects
    - j: an index to select a particular plane of data in each HDU
    
    Outputs:
    - a list of fake HDU objects
    '''
    # tqdm.write("[INFO] Creating fake HDUs for reprojecting.")
    fake_hdus = []
    for i in range(len(hdus)):
        data = hdus[i].data.copy()
        header = hdus[i].header.copy()
        data = data[j]
        del header['*3*']
        del header['*4*']
        header['WCSAXES'] = 2
        fake_hdus.append(fits.PrimaryHDU(data, header))
    return fake_hdus


def get_largest_bmaj_bmin(files):
    """
    Function to find the largest BMAJ and BMIN from a list of HDUs.

    Inputs:
    - files: 

    Outputs:
    - a tuple containing the largest BMAJ and BMIN values
    """
    
    hdu_list = [fits.open(file)[0] for file in files]

    # Initialize largest values
    largest_bmaj, largest_bmin = 0, 0

    # Loop over HDUs
    for hdu in hdu_list:
        header = hdu.header

        # Check if BMAJ and BMIN exist in the header
        if 'BMAJ' in header and 'BMIN' in header:
            bmaj = header['BMAJ']
            bmin = header['BMIN']

            # Update largest values
            if bmaj > largest_bmaj:
                largest_bmaj = bmaj

            if bmin > largest_bmin:
                largest_bmin = bmin

    return largest_bmaj, largest_bmin


#### NOT IN USE - too many files filling up /tmp/
def smooth_hdus_spectral_cube(files, largest_bmaj, largest_bmin):
    """
    Function to smooth a list of HDUs to a common resolution using SpectralCube.

    Inputs:
    - hdu_list: a list of HDU file names
    - largest_bmaj: the largest BMAJ to smooth to
    - largest_bmin: the largest BMIN to smooth to

    Outputs:
    - None, but creates smoothed images in the same directory as the input images
    """
    # Convert the largest_bmaj and largest_bmin to degrees
    largest_bmaj_deg = largest_bmaj *1.1
    largest_bmin_deg = largest_bmin *1.1

    # Create a beam object for the desired resolution
    target_beam = Beam(major=largest_bmaj_deg*u.deg, minor=largest_bmaj_deg*u.deg, pa=0*u.deg)

    hdu_list = [fits.open(file)[0] for file in files]

    # Loop over HDU file names
    for i in range(len(files)):

        hdu = hdu_list[i]

        # Load the spectral cube
        cube = SpectralCube.read(hdu)
        cube.allow_huge_operations=True

        # Convolve the cube to the target resolution
        cube = cube.convolve_to(target_beam)

        # Write the convolved cube to a new FITS file
        tqdm.write('INFO Smoothed to largest beam.')

        outfile = files[i].replace('.fits', '.smoothed.fits')
        tqdm.write('INFO Save smoothed files: %s '% outfile)
        cube.write(outfile, overwrite=True)

        del cube
        gc.collect()

    return()


def smooth_hdus_casa(files, largest_bmaj, largest_bmin):
    """
    Function to smooth a list of HDUs to a common resolution using CASA's imsmooth.

    Inputs:
    - hdu_list: a list of HDU file names
    - largest_bmaj: the largest BMAJ to smooth to (in degrees)
    - largest_bmin: the largest BMIN to smooth to (in degrees)

    Outputs:
    - None, but creates smoothed images in the same directory as the input images
    """
    # Convert the largest_bmaj and largest_bmin to degrees and increase by 10% for safety
    largest_bmaj_deg = largest_bmaj * 1.1
    largest_bmin_deg = largest_bmin * 1.1

    # Loop over HDU file names
    for file in tqdm(files, desc = 'Smoothing', position=1, leave=True):

        # Convert FITS to CASA Image
        casa_image = file.replace('.fits', '.image')
        importfits(fitsimage=file, imagename=casa_image, overwrite=True)

        # Smooth the CASA Image
        smoothed_image = casa_image + ".smoothed"
        imsmooth(imagename=casa_image, kernel='gauss', 
                major=str(largest_bmaj_deg)+'deg', minor=str(largest_bmaj_deg)+'deg', 
                pa='0deg', targetres=True, outfile=smoothed_image, overwrite=True)

        # Convert smoothed CASA Image back to FITS
        outfile = file.replace('.fits', '.smoothed.fits')
        exportfits(imagename=smoothed_image, fitsimage=outfile, overwrite=True)

        # Remove intermediate CASA images
        # os.remove(casa_image)
        # os.remove(smoothed_image)
        rmtables([casa_image, smoothed_image])

        tqdm.write('INFO Smoothed to largest beam.')
        tqdm.write('INFO Saved smoothed file: %s '% outfile)

    return


def rebin_hdus_casa(files, factor=3, overwrite=True):
    """

    """

    # Loop over HDU file names
    for file in tqdm(files, desc = 'Rebinning', position=1, leave=True):

        input_image = file.replace('.fits', '.image')
        regrid_image = file.replace('.fits', '.rebin.image')
        output_fits = file.replace('.fits', '.rebin.fits')

        # Remove any pre-existing intermediate images
        rmtables([input_image, regrid_image])

        # Import the .fits files into CASA images
        importfits(fitsimage=file, imagename=input_image, overwrite=overwrite)

        # Rebin the image by some factor 
        imrebin(imagename=input_image, outfile=regrid_image, factor=[factor,factor,1], overwrite=True)

        # Export the regridded image to a .fits file
        exportfits(imagename=regrid_image, fitsimage=output_fits, overwrite=overwrite, velocity=True)

        # Clean up by removing the intermediate images
        rmtables([input_image, regrid_image])

    return


def create_smoothed_regridded_mosaics(ACES_WORKDIR, MOLECULE):
    """
    Function to create a weighted mosaic of the TP+7m+12m cubes, if it does not exist already.

    Inputs:
    - ACES_WORKDIR: A Path object pointing to the ACES working directory.
    - MOLECULE: A string indicating the molecule for which the weighted mosaic should be created.

    Outputs:
    - None, but writes a FITS file of the weighted mosaic to the ACES working directory if it does not exist already.
    """
    # Find all the FITS files for the TP+7m+12m cubes and the 12m weights
    TP_7M_12M_cube_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.TP_7M_12M_feather_all.{MOLECULE}.image.fits')]
    # TP_7M_12M_cube_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.TP_7M_feather_all.{MOLECULE}.image.fits')]
    TWELVE_M_weight_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.12M.{MOLECULE}.image.weight.fits')]

    # If the weighted mosaic of the TP+7m+12m cubes does not exist, create it
    outputfile = ACES_WORKDIR / f'{MOLECULE}.TP_7M_12M_weighted_mosaic.fits'

    # Get the largest BMAJ and BMIN
    largest_bmaj, largest_bmin = get_largest_bmaj_bmin(TP_7M_12M_cube_files)
    smooth_hdus_casa(TP_7M_12M_cube_files, largest_bmaj, largest_bmin)

    TP_7M_12M_cube_files_smooth = [filename.replace('.fits', '.smoothed.fits') for filename in TP_7M_12M_cube_files]
    TP_7M_12M_cube_files_rebin = [filename.replace('.fits', '.smoothed.rebin.fits') for filename in TP_7M_12M_cube_files]
    TWELVE_M_weight_files_rebin = [filename.replace('.fits', '.rebin.fits') for filename in TWELVE_M_weight_files]

    # Rebin data by factor of 3
    rebin_hdus_casa(TP_7M_12M_cube_files_smooth)
    rebin_hdus_casa(TWELVE_M_weight_files)

    return()

