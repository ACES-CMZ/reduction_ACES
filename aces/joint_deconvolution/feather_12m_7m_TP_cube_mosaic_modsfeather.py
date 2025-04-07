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
from casatasks import imhead, exportfits, imtrans, feather, imreframe
from astropy.table import Table
from tqdm.auto import tqdm
import glob 


def check_files_exist(file_names):
    '''
    Function to check if all files in a list exist.
    
    Inputs:
    - file_names: a list of file paths as strings
    
    Outputs:
    - a boolean value indicating whether all files exist (True) or not (False)
    '''
    # print("[INFO] Checking if all files exist.")
    return all(file_name is not None and Path(file_name).exists() for file_name in file_names)


def get_file(filename):
    '''
    Function to find a file matching a given filename.
    
    Input:
    - filename: a string representing the name of the file to find
    
    Output:
    - the path to the last file that matches the filename, or None if no such file is found
    '''
    # print("[INFO] Getting file that matches the filename.")
    files = glob.glob(filename)
    if len(files)==0:
        print(f"[INFO] No files matching '{filename}' were found.")
        return None
    if len(files)>1:
        files.sort()
        print(f"[INFO] Too many files matching '{filename}' were found - taking highest s number.")

    return str(files[-1])


def export_fits(imagename, fitsimage, overwrite=True):
    '''
    Function to export an image to a FITS file if it doesn't already exist.
    
    Inputs:
    - imagename: a string representing the name of the image to export
    - fitsimage: a string representing the name of the FITS file to which to export the image
    
    Outputs:
    - None, but a FITS file is created or overwritten if it already exists
    '''
    # print("[INFO] Exporting image to FITS file.")
    if ((not Path(fitsimage).exists()) or (overwrite)):
        exportfits(imagename=imagename, fitsimage=fitsimage, dropdeg=True, overwrite=True)
        print(f"[INFO] Image exported to {fitsimage}.")


def process_string(input_string):
    '''
    Function to process a string - remove spaces and convert to lowercase.
    
    Input:
    - input_string: a string to be processed
    
    Output:
    - a string with all spaces removed and all characters converted to lowercase
    '''
    # print("[INFO] Processing string.")
    return input_string.replace(' ', '').lower().replace('-', '').replace('(', '').replace(')', '')


def read_table(filename):
    """
    Function to read the CSV file into an Astropy Table and process it

    Inputs:
    - filename: A string specifying the path to the CSV file

    Outputs:
    - line_spws: A dictionary with line SPWs information
    """
    table = Table.read(filename, format='csv')
    table = table[~table['Line'].mask]
    lines = table['Line'].tolist()
    line_spws = {}
    for line in lines:
        mask = table['Line'] == line
        key = process_string(line)
        line_spws[key] = {
            "mol_12m_spw": "%i" % table['12m SPW'][mask],
            "mol_7m_spw": "%i" % table['7m SPW'][mask],
            "mol_TP_spw": "%i" % table['TP SPW'][mask]
        }
    return line_spws


def feathercubes(obs_dir, obs_id, tp_cube, seven_m_cube, twelve_m_cube, twelve_m_wt, MOLECULE):
    """
    A function to process ALMA observation data.
    
    Args:
    - obs_dir (str): directory of the observation
    - obs_id (str): identifier of the observation
    - tp_cube (str): path to the total power (TP) cube
    - seven_m_cube (str): path to the 7m cube
    - twelve_m_cube (str): path to the 12m cube
    - twelve_m_wt (str): path to the 12m weight
    - MOLECULE (str): molecule under study
    
    Returns:
    None
    """
    print(f'[INFO] Processing {obs_dir}')

    # Extract the rest frequencies from the TP and 7m cubes
    tp_freq = imhead(tp_cube, mode='get', hdkey='restfreq')
    seven_m_freq = imhead(seven_m_cube, mode='get', hdkey='restfreq')

    # If the rest frequencies do not match and the reframed TP cube does not exist, reframe the TP cube to match the 7m cube
    if tp_freq['value'] != seven_m_freq['value'] and not (
            Path(tp_cube).parent / tp_cube.replace('.fits', '.reframe')
    ).is_dir():
        imreframe(
            imagename=tp_cube,
            restfreq=seven_m_freq['value'] + ' Hz',
            output=tp_cube.replace('.fits', '.reframe')
        )

    # If the transposed TP cube does not exist, transpose the TP cube
    if not (Path(tp_cube).parent / tp_cube.replace('.fits', '.imtrans')).is_dir():
        imtrans(
            imagename=tp_cube.replace('.fits', '.reframe')
            if (Path(tp_cube).parent / tp_cube.replace('.fits', '.reframe')).is_dir() else tp_cube,
            outfile=tp_cube.replace('.fits', '.imtrans'),
            order="0132"
        )

    # If the feathered TP+7m cube does not exist, feather the 7m and TP cubes together
    if not (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image').is_dir():
        feather(
            imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image'),
            highres=seven_m_cube,
            lowres=tp_cube.replace('.fits', '.imtrans')
        )

    tp_7m_cube = str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image')

    # Extract the rest frequencies from the feathered TP+7m cube and the 12m cube
    tp_7m_cube_freq = imhead(tp_7m_cube, mode='get', hdkey='restfreq')
    twelve_m_freq = imhead(twelve_m_cube, mode='get', hdkey='restfreq')

    # If the rest frequencies do not match and the reframed TP+7m cube does not exist, reframe the TP+7m cube to match the 12m cube
    if tp_7m_cube_freq['value'] != twelve_m_freq['value'] and not (Path(tp_7m_cube) / '.reframe').is_dir():
        imreframe(
            imagename=tp_7m_cube,
            restfreq=twelve_m_freq['value'] + ' Hz',
            output=tp_7m_cube + '.reframe'
        )

    # If the feathered TP+7m+12m cube does not exist, feather the 12m and TP+7m cubes together
    if not (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image').is_dir():
        feather(
            imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image'),
            highres=twelve_m_cube,
            lowres=tp_7m_cube + '.reframe' if (Path(tp_7m_cube) / '.reframe').is_dir() else tp_7m_cube
        )

    # Export the feathered cubes and the 12m weights to FITS files
    export_fits(
        imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image'),
        fitsimage=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image.fits')
    )

    export_fits(
        imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image'),
        fitsimage=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image.fits')
    )

    export_fits(
        imagename=twelve_m_wt,
        fitsimage=str(obs_dir / f'Sgr_A_st_{obs_id}.12M.{MOLECULE}.image.weight.fits')
    )

    return()


def create_feathercubes(ACES_WORKDIR, ACES_DATA, line_spws, MOLECULE):
    """
    Iteratively constructs observation cubes for each single-beam observation for the given molecule.

    Inputs:
    - ACES_WORKDIR: Path object indicating the working directory for ACES.
    - ACES_DATA: Path object indicating the data directory for ACES.
    - line_spws: Dictionary containing spectral window (SPW) information for the molecules.
    - MOLECULE: String indicating the name of the molecule for which cubes are to be generated.

    Outputs:
    - None. However, the function creates and saves observation cubes in the specified directory.
    """
    
    # Load the SB information
    sb_names = pd.read_csv('../../tables/aces_SB_uids.csv')

    # Define parts of the file naming scheme
    generic_name = '.Sgr_A_star_sci.spw'
    prefix = 'member.uid___A001_'

    # Loop over each single-beam observation (EB = Element Block)
    for i in tqdm(range(len(sb_names)), desc = 'EBs'):

        # Extract the observation ID
        obs_id = sb_names['Obs ID'][i]
        
        # Create the directory for the current observation if it doesn't exist
        obs_dir = ACES_WORKDIR / f'Sgr_A_st_{obs_id}'
        obs_dir.mkdir(exist_ok=True) 

        # Extract the MOUS IDs (Member of the Observing Unit Set) for different types of observation
        tp_mous_id = sb_names['TP MOUS ID'][i]
        seven_m_mous_id = sb_names['7m MOUS ID'][i]
        twelve_m_mous_id = sb_names['12m MOUS ID'][i]

        # Define file paths for the 7m, total power (TP), and 12m cubes and the 12m weight
        tp_cube = get_file(f"{ACES_DATA / (prefix + tp_mous_id) / 'product'}/*{generic_name}{line_spws[MOLECULE]['mol_TP_spw']}.cube.I.sd.{MOLECULE}.img")
        seven_m_cube = get_file(f"{ACES_DATA / (prefix + seven_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_7m_spw']}.cube.I.iter1.image.pbcor.{MOLECULE}.img")
        twelve_m_cube = get_file(f"{ACES_DATA / (prefix + twelve_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_12m_spw']}.cube.I.iter1.image.pbcor.{MOLECULE}.img")
        twelve_m_wt = get_file(f"{ACES_DATA / (prefix + twelve_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_12m_spw']}.cube.I.iter1.weight.{MOLECULE}.img")

        # Perform the feathering process to merge different types of observations into a single observation cube
        feathercubes(obs_dir, obs_id, tp_cube, seven_m_cube, twelve_m_cube, twelve_m_wt, MOLECULE)

    return()

        
