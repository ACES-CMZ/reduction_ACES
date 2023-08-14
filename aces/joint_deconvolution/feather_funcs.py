import glob
import pandas as pd
from pathlib import Path
from astropy.table import Table
from tqdm import tqdm
from casatasks import imhead, exportfits, imtrans, feather, imreframe, imsmooth


def check_files_exist(file_names):
    """
    Does what it says on the tin.
    """
    return all(file_name is not None and Path(file_name).exists() for file_name in file_names)


def get_file(filename):
    """
    Retrieve the file matching the provided pattern using glob.
    If multiple files are found, it sorts the files and returns the last one
    (this should not be necessary once we remove duplicate images).
    """
    files = glob.glob(filename)
    if len(files) == 0:
        print(f"[INFO] No files matching '{filename}' were found.")
        return None
    if len(files) > 1:
        files.sort()
        print(f"[INFO] Too many files matching '{filename}' were found - taking highest pipeline stage number.")

    return str(files[-1])


def export_fits(imagename, fitsimage):
    """
    Exports an image to a FITS file, dropping the degenerate axes.
    """
    if not Path(fitsimage).exists():
        exportfits(imagename=imagename, fitsimage=fitsimage, dropdeg=True)
        print(f"[INFO] Image exported to {fitsimage}.")


def process_string(input_string):
    """
    Remove spaces, dashes, and parentheses from a string and convert to lowercase.
    """
    return input_string.replace(' ', '').lower().replace('-', '').replace('(', '').replace(')', '')


def read_table(filename):
    """
    Read CSV file containing molecular line information into a dictionary.
    """
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


def feathercubes(obs_dir, obs_id, tp_cube, seven_m_cube, twelve_m_cube, twelve_m_wt, MOLECULE):
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
    twelve_m_wt : str
        Path to the 12m weight data.
    MOLECULE : str
        The molecule being processed (e.g. 'HNCO').
    """
    if check_files_exist([tp_cube, seven_m_cube, twelve_m_cube, twelve_m_wt]):
        print(f'[INFO] Processing {obs_dir}')
        if not (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image').is_dir():
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
            # This is often necessary -- feathering will not work if the axes are not in the same order
            if not (Path(tp_cube).parent / tp_cube.replace('.fits', '.imtrans')).is_dir():
                imtrans(
                    imagename=tp_cube.replace('.fits', '.reframe')
                    if (Path(tp_cube).parent / tp_cube.replace('.fits', '.reframe')).is_dir() else tp_cube,
                    outfile=tp_cube.replace('.fits', '.imtrans'),
                    order="0132"
                )

            feather_success = False
            try:
                if not (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image').is_dir():
                    feather(
                        imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image'),
                        highres=seven_m_cube,
                        lowres=tp_cube.replace('.fits', '.imtrans')
                    )
                feather_success = True
            except Exception as e:
                print(f"Feather failed with highres={seven_m_cube}: {e}")
                try:
                    imsmooth(imagename=seven_m_cube, kernel='commonbeam', targetres=True, outfile=seven_m_cube + '.commonbeam')

                    feather(
                        imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image'),
                        highres=seven_m_cube + '.commonbeam',
                        lowres=tp_cube.replace('.fits', '.imtrans')
                    )
                    feather_success = True
                except Exception as e:
                    # If the feather task failed again, print an error message and proceed
                    print(f"Feather task failed again with 7M data smoothed to a common beam. Skipping feathering: {e}")

            if feather_success:
                tp_7m_cube = str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image')
                tp_7m_cube_freq = imhead(tp_7m_cube, mode='get', hdkey='restfreq')
                twelve_m_freq = imhead(twelve_m_cube, mode='get', hdkey='restfreq')

                if tp_7m_cube_freq['value'] != twelve_m_freq['value'] and not (Path(tp_7m_cube) / '.reframe').is_dir():
                    imreframe(
                        imagename=tp_7m_cube,
                        restfreq=twelve_m_freq['value'] + ' Hz',
                        output=tp_7m_cube + '.reframe'
                    )
                # If the feathered TP+7m+12m cube does not exist, feather the 12m and TP+7m cubes together
                try:
                    feather(
                        imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image'),
                        highres=twelve_m_cube,
                        lowres=tp_7m_cube + '.reframe' if (Path(tp_7m_cube) / '.reframe').is_dir() else tp_7m_cube
                    )
                except Exception as e:
                    print(f"Feather failed with highres={twelve_m_cube}: {e}")
                    try:
                        imsmooth(imagename=twelve_m_cube, kernel='commonbeam', targetres=True, outfile=twelve_m_cube + '.commonbeam')

                        feather(
                            imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image'),
                            highres=twelve_m_cube + '.commonbeam',
                            lowres=tp_7m_cube + '.reframe' if (Path(tp_7m_cube) / '.reframe').is_dir() else tp_7m_cube
                        )
                    except Exception as e:
                        print(f"Feather task failed again with 12M data smoothed to a common beam. Skipping feathering: {e}")

        if (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image').is_dir() and not (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image.fits').is_file():
            # Export the feathered cubes and the 12m weights to FITS files
            export_fits(
                imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image'),
                fitsimage=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_12M_feather_all.{MOLECULE}.image.fits')
            )

            export_fits(
                imagename=twelve_m_wt,
                fitsimage=str(obs_dir / f'Sgr_A_st_{obs_id}.12M.{MOLECULE}.image.weight.fits')
            )
    else:
        print(f"One or more cubes do not exist for observation Sgr_A_st_{obs_id}. Skipping this one ...")


def feathercubes_without_12M(obs_dir, obs_id, tp_cube, seven_m_cube, seven_m_wt, MOLECULE):
    """
    Process TP and 7m cubes without 12m data, feather them, and export the result to FITS files.

    This function extracts rest frequencies from the TP and 7m cubes. If they do not match, the TP cube is reframed
    to match the 7m cube. The TP cube is transposed if needed. The 7m and TP cubes are then feathered together and
    exported to FITS files.

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
    seven_m_wt : str
        Path to the 7m weight data.
    MOLECULE : str
        The molecule being processed (e.g. 'HNCO').
    """
    if check_files_exist([tp_cube, seven_m_cube, seven_m_wt]):
        print(f'[INFO] Processing {obs_dir}')
        if not (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image').is_dir():
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

            try:
                feather(
                    imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image'),
                    highres=seven_m_cube,
                    lowres=tp_cube.replace('.fits', '.imtrans')
                )
            except Exception as e:
                print(f"Feather failed with highres={seven_m_cube}: {e}")
                try:
                    imsmooth(imagename=seven_m_cube, kernel='commonbeam', targetres=True, outfile=seven_m_cube + '.commonbeam')

                    feather(
                        imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image'),
                        highres=seven_m_cube + '.commonbeam',
                        lowres=tp_cube.replace('.fits', '.imtrans')
                    )
                except Exception as e:
                    print(f"Feather task failed again with 7M data smoothed to a common beam. Skipping feathering: {e}")

        if (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image').is_dir() and not (obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image.fits').is_file(): 
            export_fits(
                imagename=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image'),
                fitsimage=str(obs_dir / f'Sgr_A_st_{obs_id}.TP_7M_feather_all.{MOLECULE}.image.fits')
            )

            export_fits(
                imagename=seven_m_wt,
                fitsimage=str(obs_dir / f'Sgr_A_st_{obs_id}.7M.{MOLECULE}.image.weight.fits')
            )


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

    # Loop over each SB
    for i in tqdm(range(len(sb_names)), desc='EBs'):
        obs_id = sb_names['Obs ID'][i]
        obs_dir = ACES_WORKDIR / f'Sgr_A_st_{obs_id}'
        obs_dir.mkdir(exist_ok=True)
        tp_mous_id = sb_names['TP MOUS ID'][i]
        seven_m_mous_id = sb_names['7m MOUS ID'][i]

        seven_m_cube = get_file(
            f"{ACES_DATA / (prefix + seven_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_7m_spw']}.cube.I.iter1.image.pbcor"
        )
        tp_cube = get_file(
            f"{ACES_DATA / (prefix + tp_mous_id) / 'product'}/*{generic_name}{line_spws[MOLECULE]['mol_TP_spw']}.cube.I.sd.fits"
        )

        if process_12M:
            twelve_m_mous_id = sb_names['12m MOUS ID'][i]
            twelve_m_cube = get_file(
                f"{ACES_DATA / (prefix + twelve_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_12m_spw']}.cube.I.iter1.image.pbcor"
            )
            twelve_m_wt = get_file(
                f"{ACES_DATA / (prefix + twelve_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_12m_spw']}.cube.I.iter1.weight"
            )

            feathercubes(obs_dir, obs_id, tp_cube, seven_m_cube, twelve_m_cube, twelve_m_wt, MOLECULE)
        else:
            seven_m_wt = get_file(
                f"{ACES_DATA / (prefix + seven_m_mous_id) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_7m_spw']}.cube.I.iter1.weight"
            )

            feathercubes_without_12M(obs_dir, obs_id, tp_cube, seven_m_cube, seven_m_wt, MOLECULE)

    return ()