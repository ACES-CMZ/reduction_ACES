import os
import pandas as pd
from casatasks import feather, imregrid, imhead, immath, exportfits, importfits, imreframe


# Function to check if files exist
def check_files_exist(file_names):
    for file_name in file_names:
        if not os.path.exists(file_name):
            return False
    return True


# Get environment variables
ACES_ROOTDIR = os.getenv('ACES_ROOTDIR')
ACES_DATA = os.getenv('ACES_DATA')
ACES_WORKDIR = os.getenv('ACES_WORKDIR')

os.chdir(ACES_WORKDIR)

# Read SB uids from CSV file
sb_names = pd.read_csv(os.path.join(ACES_ROOTDIR, 'aces/data/tables/aces_SB_uids.csv'))

# Select molecule, and define generic parts of filenames
molecule = 'HNCO'
generic_name = '.Sgr_A_star_sci.spw'
prefix = 'member.uid___A001_'

# Dictionary for line spws (to be updated with other molecules and spws)
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

# Loop over all regions
for i in range(len(sb_names)):
    obs_id = sb_names['Obs ID'][i]
    obs_dir = os.path.join(ACES_WORKDIR, f'Sgr_A_st_{obs_id}')

    # Create a directory for the current region
    if not os.path.exists(obs_dir):
        os.mkdir(obs_dir)

    # Get MOUS IDs for the different arrays
    tp_mous_id = sb_names['TP MOUS ID'][i]
    seven_m_mous_id = sb_names['7m MOUS ID'][i]
    twelve_m_mous_id = sb_names['12m MOUS ID'][i]

    # Expected file paths for cubes and pb images
    tp_cube = os.path.join(ACES_DATA, prefix + tp_mous_id, 'product', prefix + tp_mous_id + generic_name + line_spws[molecule]['mol_TP_spw'] + '.cube.I.sd.fits')
    seven_m_cube = os.path.join(ACES_DATA, prefix + seven_m_mous_id, 'product', prefix + seven_m_mous_id + generic_name + line_spws[molecule]['mol_7m_spw'] + '.cube.I.pbcor.fits')
    seven_m_pb = os.path.join(ACES_DATA, prefix + seven_m_mous_id, 'product', prefix + seven_m_mous_id + generic_name + line_spws[molecule]['mol_7m_spw'] + '.cube.I.pb.fits')
    twelve_m_cube = os.path.join(ACES_DATA, prefix + twelve_m_mous_id, 'product', prefix + twelve_m_mous_id + generic_name + line_spws[molecule]['mol_12m_spw'] + '.cube.I.pbcor.fits')
    twelve_m_pb = os.path.join(ACES_DATA, prefix + twelve_m_mous_id, 'product', prefix + twelve_m_mous_id + generic_name + line_spws[molecule]['mol_12m_spw'] + '.cube.I.pb.fits')

    # Check if all files exist, and if so, feather
    if check_files_exist([tp_cube, seven_m_cube, seven_m_pb, twelve_m_cube, twelve_m_pb]) and not os.path.isdir(os.path.join(obs_dir, f'Sgr_A_st_{obs_id}.TP_7M_12M_feather.{molecule}.image')):
        seven_m_cube_image = seven_m_cube.replace('.fits', '.image')
        tp_cube_image = tp_cube.replace('.fits', '.image')
        twelve_m_cube_image = twelve_m_cube.replace('.fits', '.image')

        # Undo primary beam correction for 7m and 12m cubes
        immath(imagename=[seven_m_cube, seven_m_pb], expr='IM0*IM1', outfile=seven_m_cube_image)
        immath(imagename=[twelve_m_cube, twelve_m_pb], expr='IM0*IM1', outfile=twelve_m_cube_image)

        if not os.path.isdir(tp_cube_image):
            importfits(fitsimage=tp_cube, imagename=tp_cube_image)

        tp_freq = imhead(tp_cube_image, mode='get', hdkey='restfreq')
        seven_m_freq = imhead(seven_m_cube_image, mode='get', hdkey='restfreq')
        twelve_m_freq = imhead(twelve_m_cube_image, mode='get', hdkey='restfreq')

        # Reframe TP cube to 7m cube (if necessary)
        if tp_freq['value'] != seven_m_freq['value']:
            imreframe(imagename=tp_cube, restfreq=seven_m_freq['value'] + ' Hz', output=tp_cube.replace('.fits', '.reframe'))

        # Regrid TP cube to 7m cube
        imregrid(imagename=tp_cube.replace('.fits', '.reframe') if os.path.isdir(tp_cube + '.reframe') else tp_cube, template=seven_m_cube_image, output=tp_cube.replace('.fits', '.regrid'))

        # Regrid 12m pb cube to TP cube
        imregrid(imagename=twelve_m_pb, template=tp_cube.replace('.fits', '.regrid'), output=twelve_m_pb.replace('.fits', '.regrid'))

        # Add 12m primary beam response to TP cube
        immath(imagename=[tp_cube.replace('.fits', '.regrid'), twelve_m_pb.replace('.fits', '.regrid')], expr='IM0*IM1', outfile=tp_cube.replace('.fits', '.regrid.depb'))

        # Feather TP cube with 7m cube
        feather(imagename=os.path.join(obs_dir, f'Sgr_A_st_{obs_id}.TP_7M_feather.{molecule}.image'), highres=seven_m_cube_image, lowres=tp_cube.replace('.fits', '.regrid.depb'))

        tp_7m_cube = os.path.join(obs_dir, f'Sgr_A_st_{obs_id}.TP_7M_feather.{molecule}.image')
        tp_7m_cube_freq = imhead(tp_7m_cube, mode='get', hdkey='restfreq')

        # Reframe TP+7m cube to 12m cube (if necessary)
        if tp_7m_cube_freq['value'] != twelve_m_freq['value']:
            imreframe(imagename=tp_7m_cube_freq, restfreq=twelve_m_freq['value'] + ' Hz', output=tp_7m_cube_freq + '.reframe')

        # Regrid TP+7m cube to 12m cube
        imregrid(imagename=tp_7m_cube + '.reframe' if os.path.isdir(tp_7m_cube + '.reframe') else tp_7m_cube, template=twelve_m_cube_image, output=tp_7m_cube + '.regrid')

        # Feather TP+7m cube with 12m cube
        feather(imagename=os.path.join(obs_dir, f'Sgr_A_st_{obs_id}.TP_7M_12M_feather.{molecule}.image'), highres=twelve_m_cube_image, lowres=tp_7m_cube + '.regrid')

        # Export feathered images to FITS
        exportfits(imagename=os.path.join(obs_dir, f'Sgr_A_st_{obs_id}.TP_7M_feather.{molecule}.image'),
                   fitsimage=os.path.join(obs_dir, f'Sgr_A_st_{obs_id}.TP_7M_feather.{molecule}.image.fits'),
                   dropdeg=True)
        exportfits(imagename=os.path.join(obs_dir, f'Sgr_A_st_{obs_id}.TP_7M_12M_feather.{molecule}.image'),
                   fitsimage=os.path.join(obs_dir, f'Sgr_A_st_{obs_id}.TP_7M_12M_feather.{molecule}.image.fits'),
                   dropdeg=True)

    else:
        print(f"One or more files do not exist for observation Sgr_A_st_{obs_id}. Skipping this one ...")
