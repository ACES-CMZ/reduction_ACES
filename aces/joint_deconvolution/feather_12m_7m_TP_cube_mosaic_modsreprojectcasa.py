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
from casatasks import imhead, exportfits, imtrans, feather, imreframe, importfits, imregrid, rmtables
from spectral_cube import SpectralCube
from tqdm import tqdm 
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
    # print("[INFO] Creating fake HDUs for reprojecting.")
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
        print('INFO Smoothed to largest beam.')

        outfile = files[i].replace('.fits', '.smoothed.fits')
        print('INFO Save smoothed files: %s '% outfile)
        cube.write(outfile, overwrite=True)

        del cube
        gc.collect()

    return()


def regrid_fits_to_template(input_fits, template_fits, output_fits, overwrite=True):
    """
    A function to load a .fits file into CASA, regrid it to match a template, and save the result as another .fits file.
    
    Args:
    - input_fits (str): the path of the input .fits file
    - template_fits (str): the path of the template .fits file
    - output_fits (str): the path of the output .fits file
    - overwrite (bool): whether to overwrite existing files with the same name. Default is False.
    
    Returns:
    None
    """
    # Check if the output file already exists
    if os.path.exists(output_fits) and not overwrite:
        print("Output file already exists. Use `overwrite=True` to overwrite it.")
        return

    # Define the names of the intermediate images
    input_image = input_fits.replace('.fits', '.tmp.image')
    template_image = template_fits.replace('.fits', '.tmp.image')
    regrid_image = input_image.replace('.image', '.tmp.regrid.image')

    # Remove any pre-existing intermediate images
    rmtables([input_image, template_image, regrid_image])

    # Import the .fits files into CASA images
    importfits(fitsimage=input_fits, imagename=input_image, overwrite=overwrite)
    importfits(fitsimage=template_fits, imagename=template_image, overwrite=overwrite)

    # Regrid the image to match the template
    imregrid(imagename=input_image, template=template_image, output=regrid_image)

    # Export the regridded image to a .fits file
    exportfits(imagename=regrid_image, fitsimage=output_fits, overwrite=overwrite)

    # Clean up by removing the intermediate images
    rmtables([input_image, template_image, regrid_image])


def weighted_reproject_and_coadd(cube_files, weight_files, dir_tmp='./tmp/'):
    '''
    Function to reproject and coadd the cubes and weights.
    
    Inputs:
    - cube_files: a list of paths to the cube files
    - weight_files: a list of paths to the weight files
    
    Outputs:
    - a HDU object representing the reprojected and coadded data
    '''
    print("[INFO] Reprojecting and co-adding cubes and weights.")
    assert len(cube_files) == len(weight_files), "Mismatched number of cubes and weights."

    if not os.path.isdir(dir_tmp):
        os.mkdir(dir_tmp)

    primary_hdus = [fits.open(cube_file)[0] for cube_file in cube_files]
    weight_hdus = [fits.open(weight_file)[0] for weight_file in weight_files]

    n_hdus = len(primary_hdus)
    shape = primary_hdus[0].shape
    data_reproject = []

    fake_hdus = create_fake_hdus(primary_hdus, 0)
    wcs_out, shape_out = find_optimal_celestial_wcs(fake_hdus)
    header_out = wcs_out.to_header_string()
    hdu_out = wcs_out.to_fits()[0] 
    hdu_out.data = np.ones(shape_out)
    hdu_out.writeto('%s/hdu_out.fits' %dir_tmp, overwrite=True)

    reprojected_data, reprojected_weights = [], []

    p_bar = tqdm(range(n_hdus*2))
    p_bar.refresh()
    for i in range(n_hdus):

        print("[INFO] Processing primary_hdu[%i]" %i)

        # Load the FITS cube file
        cube = SpectralCube.read(primary_hdus[i])
        cube.allow_huge_operations=True
        cube.write('%s/cube.fits' %dir_tmp, overwrite=True)

        if i == 0: 
            keys = ['CUNIT3', 'CTYPE3', 'CRPIX3', 'CDELT3', 
            'CRVAL3', 'WCSAXES', 'SPECSYS', 'RESTFRQ',
            'BUNIT', 'BMAJ', 'BMIN', 'BPA']
            for key in keys: 
                hdu_out.header[key] = cube.header[key]

        regrid_fits_to_template('%s/cube.fits' %dir_tmp, 
                                '%s/hdu_out.fits' %dir_tmp, 
                                '%s/cube_regrid_%i.fits' %(dir_tmp, i))

        cube_regrid = SpectralCube.read('%s/cube_regrid_%i.fits' %(dir_tmp, i))
        reprojected_data.append(cube_regrid.hdu.data)

        del cube
        del cube_regrid
        gc.collect()

        p_bar.update(1)
        p_bar.refresh()

        print("[INFO] Processing weight_hdus[%i]" %i)

        # Load the FITS cube file
        cube_weight = SpectralCube.read(weight_hdus[i])
        cube_weight.allow_huge_operations=True
        cube_weight.write('%s/cube_weight.fits' %dir_tmp, overwrite=True)

        regrid_fits_to_template('%s/cube_weight.fits' %dir_tmp, 
                                '%s/hdu_out.fits' %dir_tmp, 
                                '%s/cube_weight_regrid_%i.fits' %(dir_tmp, i))

        cube_weight_regrid = SpectralCube.read('%s/cube_weight_regrid_%i.fits' %(dir_tmp, i))
        reprojected_weights.append(cube_weight_regrid.hdu.data)

        del cube_weight
        del cube_weight_regrid
        gc.collect()
        p_bar.update(1)
        p_bar.refresh()

    weighted_data = np.array(reprojected_data) * np.array(reprojected_weights)
    data_reproject = np.nansum(weighted_data, axis=0) / np.nansum(reprojected_weights, axis=0)

    hdu_reproject = fits.PrimaryHDU(data_reproject, hdu_out.header)

    return hdu_reproject


def create_weighted_mosaic(ACES_WORKDIR, MOLECULE):
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
    smooth_hdus_spectral_cube(TP_7M_12M_cube_files, largest_bmaj, largest_bmin)
    TP_7M_12M_cube_files = [filename.replace('.fits', '.smoothed.fits') for filename in TP_7M_12M_cube_files]

    print("[INFO] Creating weighted mosaic for TP+7m+12m cubes.")
    TP_7M_12M_mosaic_hdu = weighted_reproject_and_coadd(TP_7M_12M_cube_files, TWELVE_M_weight_files)
    TP_7M_12M_mosaic_hdu.writeto(outputfile, overwrite=True)
    print(f"[INFO] Created and saved weighted mosaic to {outputfile}")

    return()

