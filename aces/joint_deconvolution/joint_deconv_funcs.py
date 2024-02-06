import re
import glob
import json
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from astropy.table import Table
from casatasks import split, mstransform, tclean, importfits, exportfits, concat


def adjust_threshold(threshold_str, factor):
    match = re.match(r'([0-9.]+)([a-zA-Z]+)', threshold_str)
    if not match:
        raise ValueError("Problem with threshold parameter. Expected format: '<value><unit>'")

    numeric_part, unit = match.groups()
    new_numeric_part = float(numeric_part) / factor
    new_threshold_str = f"{new_numeric_part}{unit}"

    return new_threshold_str


def make_obs_dir(ACES_WORKDIR, obs_id):
    obs_dir = ACES_WORKDIR / f'Sgr_A_st_{obs_id}'
    obs_dir.mkdir(exist_ok=True)
    return obs_dir


def get_visfiles(ACES_DATA, prefix, mous_id, generic_name, line_spw):
    data_pattern = f"{ACES_DATA / (prefix + mous_id) / 'calibrated/working'}/*.ms"
    data_files = glob.glob(data_pattern)
    return [file for file in data_files if f"SPW{line_spw}" in file]


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


def split_visfiles(visfiles, obs_dir, LINE_SPWS, LINE, ARRAY, field='Sgr_A_star'):
    split_files = []
    for visfile in visfiles:
        outputvis = f"{obs_dir}/{Path(visfile).stem}_split.SPW{LINE_SPWS[LINE]['mol_'+ARRAY+'_spw']}.ms"
        if not Path(outputvis).exists():
            try:
                split(
                    vis=visfile,
                    spw=LINE_SPWS[LINE]['mol_'+ARRAY+'_spw'],
                    field=field,
                    outputvis=outputvis,
                    datacolumn='corrected',
                )
            except Exception as e:
                print(f"[ERROR] Splitting {visfile} failed with error: {e}. Trying again with datacolumn='data'.")
                try:
                    split(
                        vis=visfile,
                        spw=LINE_SPWS[LINE]['mol_'+ARRAY+'_spw'],
                        field=field,
                        outputvis=outputvis,
                        datacolumn='data',
                    )
                except Exception as e:
                    print(f"[ERROR] Splitting {visfile} failed with error: {e}. Something went wrong.")
            split_files.append(outputvis)
        if Path(outputvis).exists():
            split_files.append(outputvis)
    return split_files



def do_mstransform(visfiles, obs_dir, RESTFREQ, V_START, V_WIDTH, NCHAN):
    for visfile in visfiles:
        if not Path(f"{obs_dir}/{Path(visfile).stem}.mstransform").exists():
            mstransform(
                vis=visfile,
                outputvis=f"{obs_dir}/{Path(visfile).stem}.mstransform",
                outframe='LSRK',
                mode='velocity',
                veltype='radio',
                datacolumn='data',
                restfreq=str(RESTFREQ)+'GHz',
                start=V_START,
                width=V_WIDTH,
                nchan=NCHAN,
                regridms=True
            )


def do_clean(ACES_ROOTDIR, obs_dir, obs_id, tp_cube, concatvis, LINE, DEEP_CLEAN, RELAXED_MASKING, TP_STARTMODEL):

    ### THIS NEEDS TO BE UPDATED TO MERGE THE DEFAULT TCLEAN COMMANDS WITH THE OVERRIDE COMMANDS ###
    tclean_commands = json.load(open(ACES_ROOTDIR / 'aces/pipeline_scripts/default_tclean_commands.json', 'r'))
    tclean_pars = tclean_commands['Sgr_A_st_'+obs_id+'_03_TM1']['tclean_cube_pars']['spw31']
    imagename_suffix = ''

    if DEEP_CLEAN:
        tclean_pars['threshold'] = adjust_threshold(tclean_pars['threshold'], 2.0)
        imagename_suffix += '.deep_clean'

    if RELAXED_MASKING:
        tclean_pars['sidelobethreshold'] = 2.0
        tclean_pars['noisethreshold'] = 4.25
        tclean_pars['lownoisethreshold'] = 1.5
        tclean_pars['growiterations'] = 75
        imagename_suffix += '.relaxed_mask'

    """
    Commenting out for now. Automatic regridding in tclean is not working.
    Imregrid works, but requires a template image.
    Possible solution is to create dirty image, use that as template, and then regrid.
    """
    # if TP_STARTMODEL:
    #     tclean_pars['startmodel'] = tp_cube
    #     imagename_suffix += '.tp_startmodel'

    tclean_pars['vis'] = concatvis
    tclean_pars['spw'] = ''
    tclean_pars['imagename'] = str(obs_dir / f'Sgr_A_st_{obs_id}.{LINE}{imagename_suffix}')
    tclean_pars['start'] = ''
    tclean_pars['width'] = ''
    tclean_pars['nchan'] = -1
    tclean_pars['calcpsf'] = True
    tclean_pars['calcres'] = True
    tclean_pars['cyclefactor'] = 2.0
    tclean_pars['antenna'] = ''
    tclean_pars['scan'] = ''
    tclean_pars['niter'] = 0
    tclean_pars['interactive'] = False
    tclean_pars['parallel'] = True

    if not Path(tclean_pars['imagename'] + '.image').exists():
        tclean(**tclean_pars)

    if Path(tclean_pars['imagename'] + '.image.pbcor').exists() and not Path(tclean_pars['imagename'] + '.image.pbcor.fits').exists():
        export_fits(tclean_pars['imagename'] + '.image.pbcor', tclean_pars['imagename'] + '.image.pbcor.fits')

def do_joint_deconvolution(ACES_WORKDIR, ACES_DATA, ACES_ROOTDIR, REGION, LINE_SPWS, LINE, RESTFREQ, V_START, V_WIDTH, NCHAN, DEEP_CLEAN, RELAXED_MASKING, TP_STARTMODEL):
    generic_name = '.Sgr_A_star_sci.spw'
    prefix = 'member.uid___A001_'
    sb_names = pd.read_csv(ACES_ROOTDIR / 'aces/data/tables/aces_SB_uids.csv')

    for i in tqdm(range(len(sb_names)), desc='EBs'):
        obs_id = sb_names['Obs ID'][i]
        if obs_id in REGION:
            obs_dir = make_obs_dir(ACES_WORKDIR, obs_id)

            seven_m_mous_id = sb_names['7m MOUS ID'][i]
            twelve_m_mous_id = sb_names['12m MOUS ID'][i]

            tp_mous_id = sb_names['TP MOUS ID'][i]
            tp_cube = get_file(
                f"{ACES_DATA / (prefix + tp_mous_id) / 'product'}/*{generic_name}{LINE_SPWS[LINE]['mol_TP_spw']}.cube.I.sd.fits"
            )
            """
            Commenting out for now. Automatic regridding in tclean is not working.
            Imregrid works, but requires a template image.
            Possible solution is to create dirty image, use that as template, and then regrid.
            """
            # if TP_STARTMODEL:
            #     if tp_cube is None:
            #         print(f"[ERROR] No TP cube found for {obs_id}. Skipping.")
            #         continue
            #     else:
            #         importfits(fitsimage=tp_cube, imagename=tp_cube.replace('.fits', '.image'))
            #         tp_cube = tp_cube.replace('.fits', '.image')

            twelve_m_data = glob.glob(f"{ACES_DATA / (prefix + twelve_m_mous_id) / 'calibrated/working'}/*.ms")
            seven_m_data = glob.glob(f"{ACES_DATA / (prefix + seven_m_mous_id) / 'calibrated/working'}/*.ms")

            twelve_m_visfiles = split_visfiles(twelve_m_data, obs_dir, LINE_SPWS, LINE, ARRAY='12m')
            seven_m_visfiles = split_visfiles(seven_m_data, obs_dir, LINE_SPWS, LINE, ARRAY='7m')

            do_mstransform(twelve_m_visfiles, obs_dir, RESTFREQ, V_START, V_WIDTH, NCHAN)
            do_mstransform(seven_m_visfiles, obs_dir, RESTFREQ, V_START, V_WIDTH, NCHAN)

            twelve_m_visfiles_mstransform = glob.glob(str(obs_dir) + '/*split.SPW' + LINE_SPWS[LINE]['mol_12m_spw'] + '.mstransform')
            seven_m_visfiles_mstransform = glob.glob(str(obs_dir) + '/*split.SPW' + LINE_SPWS[LINE]['mol_7m_spw'] + '.mstransform')

            if (not Path(f"{obs_dir}/Sgr_A_st_{obs_id}.{LINE}.12m7m.concat.ms").exists() and
                len(seven_m_visfiles_mstransform) > 0 and
                    len(twelve_m_visfiles_mstransform) > 0):
                concat(vis=seven_m_visfiles_mstransform + twelve_m_visfiles_mstransform,
                            concatvis=f"{obs_dir}/Sgr_A_st_{obs_id}.{LINE}.12m7m.concat.ms")

            concatvis = f"{obs_dir}/Sgr_A_st_{obs_id}.{LINE}.12m7m.concat.ms"
            do_clean(ACES_ROOTDIR, obs_dir, obs_id, tp_cube, concatvis, LINE, DEEP_CLEAN, RELAXED_MASKING, TP_STARTMODEL)