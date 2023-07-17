import os
import glob
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs
from casatools import imhead, exportfits, imtrans, feather, imreframe

def check_files_exist(file_names):
    return all(file_name is not None and Path(file_name).exists() for file_name in file_names)


# Function to create fake HDUs for use in reprojecting
def create_fake_hdus(hdus, j):
    fake_hdus = []
    for i in range(len(hdus)):
        data = hdus[i].data.copy()
        header = hdus[i].header.copy()
        data = data[j]
        del header['*3*']
        del header['*4*']
        fake_hdus.append(fits.PrimaryHDU(data, header))
    return fake_hdus


# Function to reproject and coadd the cubes and weights 
def weighted_reproject_and_coadd(cube_files, weight_files):
    assert len(cube_files) == len(weight_files), "Mismatched number of cubes and weights."

    primary_hdus = [fits.open(cube_file)[0] for cube_file in cube_files]
    weight_hdus = [fits.open(weight_file)[0] for weight_file in weight_files]

    n_hdus = len(primary_hdus)
    shape = primary_hdus[0].shape
    data_reproject = []

    for j in range(shape[0]):
        fake_hdus = create_fake_hdus(primary_hdus, j)
        fake_weight_hdus = create_fake_hdus(weight_hdus, j)

        wcs_out, shape_out = find_optimal_celestial_wcs(fake_hdus)

        reprojected_data, reprojected_weights = [], []
        for i in range(n_hdus):
            array, _ = reproject_interp(fake_hdus[i], wcs_out, shape_out=shape_out)
            weight_array, _ = reproject_interp(fake_weight_hdus[i], wcs_out, shape_out=shape_out)
            reprojected_data.append(array)
            reprojected_weights.append(weight_array)

        weighted_data = np.array(reprojected_data) * np.array(reprojected_weights)
        array = np.nansum(weighted_data, axis=0) / np.nansum(reprojected_weights, axis=0)

        data_reproject.append(array)

    hdu_reproject = fits.PrimaryHDU(data_reproject, primary_hdus[0].header)

    for key in wcs_out.to_header().keys():
        if key == 'WCSAXES':
            continue
        hdu_reproject.header[key] = wcs_out.to_header()[key]

    return hdu_reproject


def get_file(filename):
    path = Path(filename)
    files = list(path.parent.glob(path.name))
    if not files:
        print(f"No files matching '{filename}' were found.")
        return None
    return str(files[-1])


def export_fits(imagename, fitsimage):
    if not Path(fitsimage).exists():
        exportfits(imagename=imagename, fitsimage=fitsimage, dropdeg=True)


# Set the relevant paths
ACES_ROOTDIR = Path(os.getenv('ACES_ROOTDIR'))
ACES_DATA = Path(os.getenv('ACES_DATA'))
ACES_WORKDIR = Path(os.getenv('ACES_WORKDIR'))
Path.cwd = ACES_WORKDIR

"""
Select the molecule/SPW based on the dictionary below, and the code will do the rest.
"""
MOLECULE = 'HNCO'

line_spws = {
    "HCOp": {
        "mol_12m_spw": "29",
        "mol_7m_spw": "20",
        "mol_TP_spw": "21"
    },
    "HNCO": {
        "mol_12m_spw": "31",
        "mol_7m_spw": "22",
        "mol_TP_spw": "23"
    },
    "cont1": {
        "mol_12m_spw": "33",
        "mol_7m_spw": "24",
        "mol_TP_spw": "25"
    },
    "cont2": {
        "mol_12m_spw": "35",
        "mol_7m_spw": "26",
        "mol_TP_spw": "27"
    }
}

sb_names = pd.read_csv(ACES_ROOTDIR / 'aces/data/tables/aces_SB_uids.csv')
generic_name = '.Sgr_A_star_sci.spw'
prefix = 'member.uid___A001_'

# Loop over all the SBs and create the 12m+7m+TP cubes for the given molecule
for i in range(len(sb_names)):
    obs_id = sb_names['Obs ID'][i]
    obs_dir = ACES_WORKDIR / f'Sgr_A_st_{obs_id}'
    obs_dir.mkdir(exist_ok=True)
    tp_mous_id = sb_names['TP MOUS ID'][i]
    seven_m_mous_id = sb_names['7m MOUS ID'][i]
    twelve_m_mous_id = sb_names['12m MOUS ID'][i]
    
    seven_m_cube = get_file(
        f"{ACES_DATA / (prefix + seven_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_7m_spw']}.cube.I.iter1.image.pbcor"
    )
    tp_cube = get_file(
        f"{ACES_DATA / (prefix + tp_mous_id) / 'product'}/*{generic_name}{line_spws[MOLECULE]['mol_TP_spw']}.cube.I.sd.fits"
    )
    twelve_m_cube = get_file(
        f"{ACES_DATA / (prefix + twelve_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_12m_spw']}.cube.I.iter1.image.pbcor"
    )
    twelve_m_wt = get_file( 
        f"{ACES_DATA / (prefix + twelve_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_12m_spw']}.cube.I.iter1.weight"
    )

    if check_files_exist([tp_cube, seven_m_cube, twelve_m_cube, twelve_m_wt]):
        print(f'Processing {obs_dir}')

        tp_freq = imhead(tp_cube, mode='get', hdkey='restfreq')
        seven_m_freq = imhead(seven_m_cube, mode='get', hdkey='restfreq')

        if tp_freq['value'] != seven_m_freq['value'] and not (
                Path(tp_cube).parent / tp_cube.replace('.fits', '.reframe')
        ).is_dir():
            imreframe(
                imagename=tp_cube,
                restfreq=seven_m_freq['value'] + ' Hz',
                output=tp_cube.replace('.fits', '.reframe')
            )

        if not (Path(tp_cube).parent / tp_cube.replace('.fits', '.imtrans')).is_dir():
            imtrans(
                imagename=tp_cube.replace('.fits', '.reframe')
                if (Path(tp_cube).parent / tp_cube.replace('.fits', '.reframe')).is_dir() else tp_cube,
                outfile=tp_cube.replace('.fits', '.imtrans'),
                order="0132"
            )

        if not (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image').is_dir():    
            feather(
                imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image'),
                highres=seven_m_cube,
                lowres=tp_cube.replace('.fits', '.imtrans')
            )

        tp_7m_cube = str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image')
        tp_7m_cube_freq = imhead(tp_7m_cube, mode='get', hdkey='restfreq')
        twelve_m_freq = imhead(twelve_m_cube, mode='get', hdkey='restfreq')

        if tp_7m_cube_freq['value'] != twelve_m_freq['value'] and not (Path(tp_7m_cube) / '.reframe').is_dir():
            imreframe(
                imagename=tp_7m_cube,
                restfreq=twelve_m_freq['value'] + ' Hz',
                output=tp_7m_cube + '.reframe'
            )

        if not (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image').is_dir():
            feather(
                imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image'),
                highres=twelve_m_cube,
                lowres=tp_7m_cube + '.reframe' if (Path(tp_7m_cube) / '.reframe').is_dir() else tp_7m_cube
            )

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

        intermediary_files = [
            tp_cube.replace('.fits', '.reframe'),
            tp_cube.replace('.fits', '.imtrans'),
            str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image'),
            tp_7m_cube + '.reframe',
            str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image')
        ]

        for file in intermediary_files:
            try:
                os.remove(file)
            except FileNotFoundError:
                print(f"File {file} not found. I cannot delete what does not exist.")

    else:
        print(f"One or more cubes do not exist for observation Sgr_A_st_{obs_id}. Skipping this one ...")

# Create a weighted cube mosaic of the TP+7m+12m data
TP_7M_12M_cube_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.TP_7M_12M_feather_all.{MOLECULE}.image.fits')]
TWELVE_M_weight_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.12M.{MOLECULE}.image.weight.fits')]

if not (ACES_WORKDIR / f'{MOLECULE}.TP_7M_12M_weighted_mosaic.fits').exists():
    TP_7M_12M_mosaic_hdu = weighted_reproject_and_coadd(TP_7M_12M_cube_files, TWELVE_M_weight_files)
    TP_7M_12M_mosaic_hdu.writeto(f'{MOLECULE}.TP_7M_12M_weighted_mosaic.fits')