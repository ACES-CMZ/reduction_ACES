# Import required modules
import os
import glob
from tqdm import tqdm
from casatasks import feather, exportfits, imtrans, imreframe, imhead, importfits
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
import numpy as np 

# Function to find files in a directory with a specific pattern
def find_files(directory, pattern):
    # Use a list comprehension to iterate over all the files in the directory and its sub-directories
    # The glob function with the '**' operator is used to search recursively in the directory
    return [f for f in glob.glob(os.path.join(directory, '**', pattern), recursive=True)]

# Function to convert the data in an HDU (Header Data Unit) to float32 format
def convert_to_float32(hdu):
    hdu.data = hdu.data.astype('float32')
    return(hdu)

# Function to crop a FITS cube to a specific velocity range
def crop_cube_velocity_range(fits_file, rest_frequency, v_start, v_end, v_res=None):

    # Load the FITS cube file
    cube = SpectralCube.read(fits_file)
    cube.allow_huge_operations=True

    # Convert the cube to velocity space using the given rest frequency
    cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=rest_frequency*u.GHz)

    # Define the velocity range for cropping
    vrange = [v_start*u.km/u.s, v_end*u.km/u.s]

    # Crop the cube without regridding
    cropped_cube = cube.spectral_slab(*vrange)

    # Out of velocity range
    if cropped_cube.shape[0]<=1: 
        return(None)

    if v_res != None: 
        # Crop the cube with regridding        
        cropped_cube = regrid_cube(cropped_cube, v_start, v_end, v_res)
    
    # Crop the cube further to the minimal enclosing subcube
    cropped_cube = cropped_cube.minimal_subcube()

    # Convert the cube to Kelvin units
    cropped_cube = cropped_cube.to(u.K)

    # Convert the cube to an HDU
    hdu = cropped_cube.hdu

    # Convert the HDU data to float32 format
    hdu = fits.PrimaryHDU(hdu.data, hdu.header)
    hdu = convert_to_float32(hdu) 

    return hdu

# Function to regrid a cube to a new velocity axis
def regrid_cube(cube, v_start, v_end, v_res):

    # Define the new velocity axis
    new_velocity = np.arange(v_start, v_end+v_res, v_res)*u.km/u.s

    # Regrid the cube to the new velocity axis
    new_cube = cube.spectral_interpolate(new_velocity, suppress_smooth_warning=True)

    return new_cube

# Function to process a string - remove spaces and convert to lowercase
def process_string(input_string):
    return input_string.replace(' ', '').lower().replace('-', '').replace('(', '').replace(')', '')



# Define path to files
path = '../rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/'

# Read the CSV file into an Astropy Table
filename = '../../tables/ACESSpectralLineCoverage.csv'
table = Table.read(filename, format='csv')

# Only keep rows in the table where the 'Line' column is not masked
table = table[~table['Line'].mask]

# Convert the 'Line' and 'Rest  (GHz)' columns to lists
lines = table['Line'].tolist()   
frequencies = table['Rest  (GHz)'].quantity.value

# Find all '.pbcor' files in the specified directory
# imagenames = find_files(path, '*.pbcor')
imagenames = find_files(path, '*.weight')
imagenames.sort()


# Loop over each image file
for imagename in tqdm(imagenames, desc='Processing images'):
    
    # Print progress information
    print('[INFO] Exporting CASA fits:', imagename.split('/')[-1])

    # Define the output FITS image filename
    fitsimage=f'{imagename}.fits'

    exportfits(imagename=imagename, 
               fitsimage=fitsimage,
               overwrite=True,
               dropstokes=True,
               history=False,
               dropdeg=True)

    # Loop over each line in the 'Line' list
    for i in range(len(lines)):

        # Define the line and its corresponding frequency
        line = lines[i]
        frequency = frequencies[i]

        # Crop the FITS cube to the specified velocity range
        hdu = crop_cube_velocity_range(fitsimage, frequency, v_start=-10, v_end=160, v_res=3)

        if hdu == None: 
            continue

        # Define the output filename
        outputfile = fitsimage.replace('.fits', '.%s.fits' %process_string(line))

        # Print progress information
        print('[INFO] Saving file:', outputfile)

        # Write the HDU to a new FITS file
        hdu.writeto(outputfile, overwrite=True)

        importfits(imagename=outputfile.replace('.fits','.img'), 
           fitsimage=outputfile,
           zeroblanks=True,
           overwrite=True)