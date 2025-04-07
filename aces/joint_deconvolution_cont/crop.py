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

def get_croppeddata(hdu, region, square=False):
    
    hdu_crop = hdu.copy()  # Create a copy of the input HDU object
    # del hdu  # Delete the original HDU object to free up memory
    # _ = gc.collect()  # Perform garbage collection

    wcs = WCS(hdu_crop)  # Create a WCS object from the HDU header

    try: 
        if square:
            radius = max([region['width'], region['height']]) # Calculate the radius for the square cutout 
            cutout = Cutout2D(hdu_crop.data, region['position'], radius, wcs=wcs)  # Create a square cutout
        else:
            cutout = Cutout2D(hdu_crop.data, region['position'], [region['width'], region['height']], wcs=wcs)  # Create a rectangular cutout
    except: 
        print('[INFO] NO cross-over')
        return(None)

    hdu_crop.data = cutout.data  # Update the data in the cropped HDU object
    hdu_crop.header.update(cutout.wcs.to_header())  # Update the header of the cropped HDU object with the cutout's WCS information

    header = hdu.header
    header_out = hdu_crop.header
    header_keywords = np.array(header.cards)[:,0]
    keys = ['BMAJ', 'BMIN', 'BPA', 'BUNIT']
    for key in keys: 
        if key in header_keywords:
            header_out[key] = header[key]

    if 'BUNIT' not in header_keywords: 
        header_out['BUNIT'] = 'Jy/beam'

    hdu_crop.header = header_out

    del cutout  # Delete the cutout to free up memory
    _ = gc.collect()  # Perform garbage collection

    return hdu_crop  # Return the cropped HDU object

def get_region(l, b, w, h, frame='galactic'): 
    region = {'position': SkyCoord(l=l *u.deg, b=b *u.deg, frame=frame), 
          'width': w *u.deg,
          'height': h *u.deg}
    return(region)

# Get a list of all .fits files in the 'data' directory
files = ['../data/processed/cont_tp.fits', '../data/processed/cont_12m.fits']

# Define the galactic coordinates and the dimensions of the desired region in the sky
l = 0.4697169  # galactic longitude
b = -0.0125704  # galactic latitude
width = 1.1808514
height = 1.1487904

# Get the defined region 
region = get_region(l, b, height, width)

# Loop over each file. We're using tqdm to create a progress bar. 
# tqdm automatically determines the total number of iterations from the length of the 'files' list.
for file in tqdm(files, desc="Processing files", unit="file"):
    
    # Print the name of the current file being processed
    print('[INFO] infile - %s' % file)

    # Open the FITS file
    hdu = fits.open(file)[0]

    # Crop the data in the FITS file to the defined region 
    hdu_crop = get_croppeddata(hdu, region)

    # Create a WCS (World Coordinate System) object from the cropped HDU
    wcs = WCS(hdu_crop)

    # Calculate the middle pixel coordinates of the cropped image
    y_mid, x_mid = hdu_crop.shape[0]/2, hdu_crop.shape[1]/2

    # Convert the middle pixel coordinates to galactic coordinates
    l_mid = wcs.pixel_to_world(x_mid, y_mid).l
    b_mid = wcs.pixel_to_world(x_mid, y_mid).b

    # Update the header information of the cropped HDU
    hdu_crop.header['BUNIT'] = 'Jy/beam'  # unit of the data
    hdu_crop.header['CRPIX1'] = x_mid  # x-coordinate of the reference pixel
    hdu_crop.header['CRPIX2'] = y_mid  # y-coordinate of the reference pixel
    hdu_crop.header['CRVAL1'] = l_mid.value  # galactic longitude of the reference pixel
    hdu_crop.header['CRVAL2'] = b_mid.value  # galactic latitude of the reference pixel

    # If the cropping was successful, save the cropped data to a new FITS file
    if hdu_crop is not None: 
        # Construct the name of the output file by replacing '.fits' with '_cropped.fits' in the input file name
        output_file = file.replace('.fits', '_cropped.fits')
        # Write the cropped HDU to the output file, overwriting it if it already exists
        hdu_crop.writeto(output_file, overwrite=True)
