import os
import glob
import shutil
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from astropy.table import Table
from casatasks import importfits, exportfits, imhead, imreframe, imtrans, feather, imsmooth


def check_files_exist(file_names):
    return all(file_name is not None and Path(file_name).exists() for file_name in file_names)


def get_file(filename):
    """
    Retrieve the file matching the provided pattern using glob.
    If multiple files are found, it sorts the files and returns the first one
    (this should not be necessary once we remove duplicate images).
    """
    files = glob.glob(filename)
    if len(files) == 0:
        print(f"[INFO] No files matching '{filename}' were found.")
        return None
    if len(files) > 1:
        files.sort()
        print(f"[INFO] Too many files matching '{filename}' were found - taking lowest pipeline stage number.")

    return str(files[0])


def export_fits(imagename, fitsimage):
    if not Path(fitsimage).exists():
        exportfits(imagename=imagename, fitsimage=fitsimage, dropdeg=True)
        print(f"[INFO] Image exported to {fitsimage}.")


def import_fits(fitsimage, imagename, overwrite=False):
    if not Path(imagename).exists() or overwrite:
        importfits(fitsimage=fitsimage, imagename=imagename, defaultaxes=True, defaultaxesvalues=['', '', '', 'I'], overwrite=overwrite)
        print(f"[INFO] FITS file imported into {imagename}.")

    return imagename


def process_string(input_string):
    return input_string.replace(' ', '').lower().replace('-', '').replace('(', '').replace(')', '')


def read_table(filename):
    table = Table.read(filename, format='csv')
    lines = table['Line'].tolist()
    line_spws = {}
    for line in lines:
        mask = table['Line'] == line
        key = process_string(line)
        line_spws[key] = {
            "mol_12m_spw": "%i" % table['12m SPW'][mask][0],
            "mol_7m_spw": "%i" % table['7m SPW'][mask][0],
            "mol_TP_spw": "%i" % table['TP SPW'][mask][0],
            "restfreq": table['Rest (GHz)'][mask][0],
        }
    return line_spws


def feathercubes(obs_dir, obs_id, tp_cube, seven_m_cube, twelve_m_cube, MOLECULE):
    """
    Process TP and 7m cubes with 12m data, feather them, and export the result to FITS files.

    This function also checks the rest frequencies of the cubes. If they do not match, the cube is reframed
    to match the target cube. The TP cube is transposed if needed.

    Occasionally, feathering will fail if one of the cubes has multiple beams.
    In this case, the offending cube is smoothed to a common beam and feathering is attempted again.

    Parameters
    ----------
    obs_dir : pathlib.Path
        The directory for the relevant region (e.g. 'Sgr_A_st_a').
    obs_id : str
        The ID of the observation (e.g. 'a').
    tp_cube : str
        Path to the total power data cube.
    seven_m_cube : str
        Path to the 7m data cube.
    twelve_m_cube : str
        Path to the 12m data cube.
    MOLECULE : str
        The molecule being processed (e.g. 'HNCO').
    """
    if check_files_exist([tp_cube, seven_m_cube, twelve_m_cube]):
        print(f'[INFO] Processing {obs_dir}')
        OUTPUT_TP_7M_12M = obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image.statcont.contsub'
        OUTPUT_TP_7M = obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image.statcont.contsub'
        if not OUTPUT_TP_7M_12M.is_dir():
            # Import and set observatory for 7m data
            seven_m_cube_im = import_fits(seven_m_cube, seven_m_cube.replace('.fits', '.image'))
            imhead(imagename=seven_m_cube_im, mode='put', hdkey='telescope', hdvalue='ALMA')
            seven_m_head = imhead(seven_m_cube_im)
            if 'perplanebeams' in seven_m_head:
                if not os.path.isdir(seven_m_cube_im + '.commonbeam'):
                    imsmooth(imagename=seven_m_cube_im, kernel='commonbeam', targetres=True, outfile=seven_m_cube_im + '.commonbeam')

            # Import and set observatory for 12m data
            twelve_m_cube_im = import_fits(twelve_m_cube, twelve_m_cube.replace('.fits', '.image'), overwrite=True)
            imhead(imagename=twelve_m_cube_im, mode='put', hdkey='telescope', hdvalue='ALMA')
            twelve_m_head = imhead(twelve_m_cube_im)
            if 'perplanebeams' in twelve_m_head:
                if not os.path.isdir(twelve_m_cube_im + '.commonbeam'):
                    imsmooth(imagename=twelve_m_cube_im, kernel='commonbeam', targetres=True, outfile=twelve_m_cube_im + '.commonbeam')

            # Set observatory for TP cube and get frequencies
            imhead(imagename=tp_cube, mode='put', hdkey='telescope', hdvalue='ALMA')
            tp_freq = imhead(tp_cube, mode='get', hdkey='restfreq')
            seven_m_freq = imhead(seven_m_cube_im, mode='get', hdkey='restfreq')

            # If the rest frequencies do not match and the reframed TP cube does not exist, reframe the TP cube to match the 7m cube
            if tp_freq['value'] != seven_m_freq['value']:
                if os.path.isdir(tp_cube.replace('.fits', '.reframe')):
                    shutil.rmtree(tp_cube.replace('.fits', '.reframe'))
                imreframe(
                    imagename=tp_cube,
                    restfreq=seven_m_freq['value'] + ' Hz',
                    output=tp_cube.replace('.fits', '.reframe')
                )

            # If the transposed TP cube does not exist, transpose the TP cube
            # This is often necessary -- feathering will not work if the axes are not in the same order
            if os.path.isdir(tp_cube.replace('.fits', '.imtrans')):
                shutil.rmtree(tp_cube.replace('.fits', '.imtrans'))
            imtrans(
                imagename=tp_cube.replace('.fits', '.reframe')
                if os.path.isdir(tp_cube.replace('.fits', '.reframe')) else tp_cube,
                outfile=tp_cube.replace('.fits', '.imtrans'),
                order="0132"
            )

            if not os.path.isdir(OUTPUT_TP_7M):
                feather_success = False
                try:
                    feather(
                        imagename=str(OUTPUT_TP_7M),
                        highres=seven_m_cube_im + '.commonbeam' if os.path.isdir(seven_m_cube_im + '.commonbeam') else seven_m_cube_im,
                        lowres=tp_cube.replace('.fits', '.imtrans')
                    )
                    feather_success = True
                except Exception as e:
                    print(f"Feather failed with highres={seven_m_cube_im}: {e}")
            else:
                feather_success = True

            if feather_success:
                tp_7m_cube = str(OUTPUT_TP_7M)
                tp_7m_cube_freq = imhead(tp_7m_cube, mode='get', hdkey='restfreq')
                twelve_m_freq = imhead(twelve_m_cube_im, mode='get', hdkey='restfreq')

                if tp_7m_cube_freq['value'] != twelve_m_freq['value']:
                    if os.path.isdir(tp_7m_cube + '.reframe'):
                        shutil.rmtree(tp_7m_cube + '.reframe')
                    imreframe(
                        imagename=tp_7m_cube,
                        restfreq=twelve_m_freq['value'] + ' Hz',
                        output=tp_7m_cube + '.reframe'
                    )
                # If the feathered TP+7m+12m cube does not exist, feather the 12m and TP+7m cubes together
                try:
                    feather(
                        imagename=str(OUTPUT_TP_7M_12M),
                        highres=twelve_m_cube_im + '.commonbeam' if os.path.isdir(twelve_m_cube_im + '.commonbeam') else twelve_m_cube_im,
                        lowres=tp_7m_cube + '.reframe' if (Path(tp_7m_cube) / '.reframe').is_dir() else tp_7m_cube
                    )
                except Exception as e:
                    print(f"Feather failed with highres={twelve_m_cube_im}: {e}")
                    try:
                        feather(
                            imagename=str(OUTPUT_TP_7M_12M),
                            highres=twelve_m_cube_im + '.commonbeam',
                            lowres=tp_7m_cube + '.reframe' if (Path(tp_7m_cube) / '.reframe').is_dir() else tp_7m_cube
                        )
                    except Exception as e:
                        print(f"Feather task failed again with 12M data smoothed to a common beam. Skipping feathering: {e}")

        if OUTPUT_TP_7M_12M.is_dir() and not os.path.isfile(str(OUTPUT_TP_7M_12M) + '.fits'):
            # Export the feathered cubes and the 12m weights to FITS files
            export_fits(
                imagename=str(OUTPUT_TP_7M_12M),
                fitsimage=str(OUTPUT_TP_7M_12M) + '.fits'
            )
    else:
        print(f"One or more cubes do not exist for observation Sgr_A_st_{obs_id}. Skipping this one ...")


def create_feathercubes(ACES_WORKDIR, ACES_DATA, ACES_ROOTDIR, line_spws, MOLECULE, process_12M=True):
    """
    Loop over each SB and identify the relevant MOUS IDs
    Retrieves the corresponding data cubes and feathers them

    The function can optionally include 12m data

    Parameters
    ----------
    ACES_WORKDIR : Path
        Path object representing the working directory.
    ACES_DATA : Path
        Path object representing the directory where the data files are stored.
    ACES_ROOTDIR : Path
        Path object representing the root directory of ACES.
    line_spws : dict
        Dictionary mapping molecule identifiers to spectral window information.
    MOLECULE : str
        The molecule being processed (e.g. 'HNCO').
    process_12M : bool, optional
        If True, 12m data is processed and included in the feathering process. Default is True.
    """
    # Load the SB information
    sb_names = pd.read_csv(ACES_ROOTDIR / 'aces/data/tables/aces_SB_uids.csv')

    # Define common parts of the file naming scheme
    generic_name = '.Sgr_A_star_sci.spw'
    prefix = 'member.uid___A001_'
    regions = ['ap']

    # Loop over each SB
    for i in tqdm(range(len(sb_names)), desc='EBs'):
        obs_id = sb_names['Obs ID'][i]
        if obs_id in regions:
            obs_dir = ACES_WORKDIR / f'Sgr_A_st_{obs_id}'
            obs_dir.mkdir(exist_ok=True)
            tp_mous_id = sb_names['TP MOUS ID'][i]
            seven_m_mous_id = sb_names['7m MOUS ID'][i]
            TM_SPW = 'SPW_' + line_spws[MOLECULE]['mol_12m_spw']

            seven_m_cube = get_file(
                f"{ACES_DATA / (prefix + seven_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_7m_spw']}.cube.I.iter1.image.pbcor.statcont.contsub.fits"
            )
            tp_cube = get_file(
                f"{ACES_DATA / (prefix + tp_mous_id) / 'product'}/*{generic_name}{line_spws[MOLECULE]['mol_TP_spw']}.cube.I.sd.fits"
            )

            if process_12M:
                twelve_m_mous_id = sb_names['12m MOUS ID'][i]
                twelve_m_cube = get_file(
                    f"{ACES_DATA / (prefix + twelve_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_12m_spw']}.cube.I.iter1.image.pbcor.statcont.contsub.fits"
                )

                feathercubes(obs_dir, obs_id, tp_cube, seven_m_cube, twelve_m_cube, TM_SPW)


###_______________________________________________________________________________________________________________________
ACES_ROOTDIR = Path(os.getenv('ACES_ROOTDIR'))
ACES_WORKDIR = Path(os.getenv('ACES_WORKDIR'))
ACES_DATA = Path(os.getenv('ACES_DATA'))

LINE_TABLE = (ACES_ROOTDIR / 'aces/data/tables/linelist.csv')
LINES = ['cs21']
INCLUDE_12M = True

for LINE in tqdm(LINES, desc='LINES'):
    LINE_SPWS = read_table(LINE_TABLE)
    create_feathercubes(ACES_WORKDIR, ACES_DATA, ACES_ROOTDIR, LINE_SPWS, LINE, process_12M=INCLUDE_12M)