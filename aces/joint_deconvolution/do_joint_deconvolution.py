import os
import re
import glob
import json
import shutil
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from astropy.table import Table
from casatasks import split, mstransform, tclean, feather, exportfits, immath, importfits, imregrid, imtrans, ia

ACES_ROOTDIR = None
ACES_WORKDIR = None
ACES_DATA = None
LINE_TABLE = None
LINE_SPWS = None


def init(aces_rootdir, aces_workdir, aces_data, line_table):
    global ACES_ROOTDIR, ACES_WORKDIR, ACES_DATA, LINE_TABLE, LINE_SPWS
    ACES_ROOTDIR = Path(aces_rootdir)
    ACES_WORKDIR = Path(aces_workdir)
    ACES_DATA = Path(aces_data)
    LINE_TABLE = ACES_ROOTDIR / line_table
    LINE_SPWS = read_table(LINE_TABLE)


def adjust_threshold(threshold_str, factor):
    match = re.match(r'([0-9.]+)([a-zA-Z]+)', threshold_str)
    if not match:
        raise ValueError("Problem with threshold parameter. Expected format: '<value><unit>'")
    numeric_part, unit = match.groups()
    new_numeric_part = float(numeric_part) / factor
    return f"{new_numeric_part}{unit}"


def make_obs_dir(obs_id):
    obs_dir = ACES_WORKDIR / f'Sgr_A_st_{obs_id}'
    obs_dir.mkdir(exist_ok=True)
    return obs_dir


def get_file(filename):
    files = glob.glob(filename)
    if not files:
        print(f"[INFO] No files matching '{filename}' were found.")
        return None
    if len(files) > 1:
        files.sort()
        print(f"[INFO] Multiple files matching '{filename}' were found - taking highest pipeline stage number.")
    return str(files[-1])


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


def split_visfiles(visfiles, obs_dir, line, array, field='Sgr_A_star'):
    split_files = []
    for visfile in visfiles:
        outputvis = f"{obs_dir}/{Path(visfile).stem}_split.SPW{LINE_SPWS[line]['mol_' + array + '_spw']}.ms"
        if not Path(outputvis).exists():
            try:
                split(
                    vis=visfile,
                    spw=LINE_SPWS[line]['mol_' + array + '_spw'],
                    field=field,
                    outputvis=outputvis,
                    datacolumn='corrected',
                )
            except Exception as e:
                print(f"[ERROR] Splitting {visfile} failed with error: {e}. Trying again with datacolumn='data'.")
                try:
                    split(
                        vis=visfile,
                        spw=LINE_SPWS[line]['mol_' + array + '_spw'],
                        field=field,
                        outputvis=outputvis,
                        datacolumn='data',
                    )
                except Exception as e:
                    print(f"[ERROR] Splitting {visfile} failed with error: {e}. Skipping this file.")
                    continue
        if Path(outputvis).exists():
            split_files.append(outputvis)
    return split_files


def do_mstransform(visfiles, obs_dir, restfreq, v_start, v_width, nchan):
    for visfile in visfiles:
        if not Path(f"{obs_dir}/{Path(visfile).stem}.mstransform").exists():
            mstransform(
                vis=visfile,
                outputvis=f"{obs_dir}/{Path(visfile).stem}.mstransform",
                outframe='LSRK',
                mode='velocity',
                veltype='radio',
                datacolumn='data',
                restfreq=f"{restfreq}GHz",
                start=v_start,
                width=v_width,
                nchan=nchan,
                regridms=True
            )


def create_dirty_image(obs_dir, obs_id, concatvis, line, tclean_pars):
    dirty_imagename = str(obs_dir / f'Sgr_A_st_{obs_id}.{line}.dirty')
    if not Path(dirty_imagename + '.image').exists():
        print(f"Creating dirty image: {dirty_imagename}")
        dirty_pars = tclean_pars.copy()
        dirty_pars.update({
            'vis': concatvis,
            'imagename': dirty_imagename,
            'spw': '',
            'start': '',
            'width': '',
            'nchan': -1,
            'antenna': '',
            'scan': '',
            'niter': 0,
            'interactive': False,
            'calcpsf': True,
            'calcres': True,
            'pbcor': False,
            'parallel': False
        })
        tclean(**dirty_pars)
    return dirty_imagename + '.image'


def convert_jy_beam_to_jy_pixel(input_image, output_image):
    if not Path(output_image).exists():
        print(f"Converting image from Jy/beam to Jy/pixel: {output_image}")
        try:
            ia.open(input_image)
            beam_info = ia.restoringbeam()
            if not beam_info:
                raise ValueError("No beam information found in the image.")
            
            bmaj = beam_info['major']['value']
            bmin = beam_info['minor']['value']
            
            if beam_info['major']['unit'] == 'arcsec':
                pass
            elif beam_info['major']['unit'] == 'deg':
                bmaj *= 3600
                bmin *= 3600
            else:
                raise ValueError(f"Unexpected beam unit: {beam_info['major']['unit']}")
            
            csys = ia.coordsys()
            increment = csys.increment()['numeric']
            pixel_area = abs(increment[0] * increment[1]) * (3600 ** 2) * ((180/np.pi) ** 2)
            ia.close()
            
            beam_area = (np.pi * bmaj * bmin) / (4 * np.log(2))
            conversion_factor = pixel_area / beam_area
            
            immath(imagename=input_image,
                   mode='evalexpr',
                   expr=f'iif(IM0 > 0, IM0 * {conversion_factor}, 0)',
                   outfile=output_image)
            
            print(f"Converted image saved as: {output_image}")
        except Exception as e:
            print(f"Error converting image to Jy/pixel: {str(e)}")
            raise

    return output_image


def apply_pb_attenuation(tp_image, pb_image, output_image):
    if not Path(output_image).exists():
        print(f"Applying primary beam attenuation: {output_image}")
        try:
            immath(imagename=[tp_image, pb_image],
                   mode='evalexpr',
                   expr='IM0*IM1',
                   outfile=output_image)
            print(f"PB-attenuated image saved as: {output_image}")
        except Exception as e:
            print(f"Error applying PB attenuation: {str(e)}")
            raise

    return output_image


def regrid_and_reorder_cube(input_image, template_image, output_image):
    if not Path(output_image).exists():
        print(f"Regridding & reordering cube: {output_image}")

        try:
            if input_image.endswith('.fits'):
                importfits(fitsimage=input_image, imagename=input_image + '.temp.image', overwrite=True)
                input_image = input_image + '.temp.image'

            regridded_image = input_image + '.regrid.image'
            imregrid(
                imagename=input_image,
                template=template_image,
                output=regridded_image,
                overwrite=True
            )

            imtrans(
                imagename=regridded_image,
                outfile=output_image,
                order='0132'
            )

            print(f"Regridded and reordered cube saved as: {output_image}")

            for temp_image in [input_image + '.temp.image', regridded_image]:
                if os.path.exists(temp_image):
                    shutil.rmtree(temp_image)

        except Exception as e:
            print(f"Error in regrid_and_reorder_cube: {str(e)}")
            raise

    return output_image


# TODO: Expand this to handle other files (e.g. mstransformed files, modified TP cubes, etc.)
def remove_existing_clean_outputs(obs_dir, obs_id, line, imagename_suffix):
    pattern = f'Sgr_A_st_{obs_id}.{line}{imagename_suffix}*'
    for item in obs_dir.glob(pattern):
        if item.is_dir():
            shutil.rmtree(item)
        else:
            item.unlink()
    print(f"Removed existing clean outputs matching: {pattern}")


# TODO: Add some of the parameters to main function to avoid hardcoding
def do_clean(obs_dir, obs_id, tp_cube, concatvis, line, deep_clean, relaxed_masking, tp_startmodel, overwrite_clean):
    # TODO: Add a function to merge the default tclean commands with the override parameters
    tclean_commands = json.load(open(ACES_ROOTDIR / 'aces/pipeline_scripts/default_tclean_commands.json', 'r'))
    tclean_pars = tclean_commands[f'Sgr_A_st_{obs_id}_03_TM1']['tclean_cube_pars']['spw31']
    imagename_suffix = ''

    if deep_clean:
        tclean_pars['threshold'] = adjust_threshold(tclean_pars['threshold'], 1.50)
        imagename_suffix += '.deep_clean'

    if relaxed_masking:
        tclean_pars['sidelobethreshold'] = 2.0
        tclean_pars['noisethreshold'] = 4.25
        tclean_pars['lownoisethreshold'] = 1.5
        tclean_pars['growiterations'] = 100
        tclean_pars['pblimit'] = 0.3
        imagename_suffix += '.relaxed_mask'

    dirty_image = create_dirty_image(obs_dir, obs_id, concatvis, line, tclean_pars)

    reordered_tp_cube = regrid_and_reorder_cube(
        tp_cube, 
        dirty_image, 
        str(obs_dir / f'tp_cube_{obs_id}.{line}.regridded.image')
    )

    if tp_startmodel:
        tp_cube_jy_per_pixel = convert_jy_beam_to_jy_pixel(
            reordered_tp_cube,
            str(obs_dir / f'tp_cube_{obs_id}.{line}.jy_per_pixel.image')
        )

        pb_image = str(obs_dir / f'Sgr_A_st_{obs_id}.{line}.dirty.pb')
        tp_cube_pb_attenuated = apply_pb_attenuation(
            tp_cube_jy_per_pixel,
            pb_image,
            str(obs_dir / f'tp_cube_{obs_id}.{line}.jy_per_pixel.pb_attenuated.image')
        )

        tclean_pars['startmodel'] = tp_cube_pb_attenuated
        imagename_suffix += '.tp_startmodel'

    tclean_pars.update({
        'vis': concatvis,
        'spw': '',
        'imagename': str(obs_dir / f'Sgr_A_st_{obs_id}.{line}{imagename_suffix}'),
        'start': '',
        'width': '',
        'nchan': -1,
        'calcpsf': True,
        'calcres': True,
        'cyclefactor': 3.0,
        'deconvolver': 'multiscale',
        'scales': [0, 8, 16, 32, 64],
        'antenna': '',
        'scan': '',
        'niter': 1000000,
        'interactive': False,
        'parallel': False,
        'pbcor': True
    })

    if overwrite_clean:
        remove_existing_clean_outputs(obs_dir, obs_id, line, imagename_suffix)

    if not Path(tclean_pars['imagename'] + '.image').exists() or overwrite_clean:
        print(f"Running tclean for {tclean_pars['imagename']}")
        tclean(**tclean_pars)

    # TODO: abstract feathering to a separate function
    highres_image = tclean_pars['imagename'] + '.image.pbcor'
    lowres_image = reordered_tp_cube
    feathered_image = tclean_pars['imagename'] + '.image.pbcor.feather'

    if not Path(feathered_image).exists() or overwrite_clean:
        try:
            feather(imagename=feathered_image, highres=highres_image, lowres=lowres_image)
            print(f"Feathering completed: {feathered_image}")
        except Exception as e:
            print(f"Feathering failed: {str(e)}")

    # TODO: abstract FITS export to a separate function
    for image_type in ['.image.pbcor', '.image.pbcor.feather']:
        input_image = tclean_pars['imagename'] + image_type
        output_fits = input_image + '.fits'
        if (Path(input_image).exists() and not Path(output_fits).exists()) or overwrite_clean:
            try:
                exportfits(imagename=input_image, fitsimage=output_fits, dropdeg=True)
                print(f"Exported FITS: {output_fits}")
            except Exception as e:
                print(f"Failed to export FITS for {input_image}: {str(e)}")


def do_joint_deconvolution(region, line, restfreq, v_start, v_width, nchan, deep_clean, relaxed_masking, tp_startmodel, overwrite_clean):
    generic_name = '.Sgr_A_star_sci.spw'
    prefix = 'member.uid___A001_'
    sb_names = pd.read_csv(ACES_ROOTDIR / 'aces/data/tables/aces_SB_uids.csv')

    for i in tqdm(range(len(sb_names)), desc='EBs'):
        obs_id = sb_names['Obs ID'][i]
        if obs_id not in region:
            continue

        obs_dir = make_obs_dir(obs_id)

        seven_m_mous_id = sb_names['7m MOUS ID'][i]
        twelve_m_mous_id = sb_names['12m MOUS ID'][i]

        tp_mous_id = sb_names['TP MOUS ID'][i]
        tp_cube = get_file(
            f"{ACES_DATA / (prefix + tp_mous_id) / 'product'}/*{generic_name}{LINE_SPWS[line]['mol_TP_spw']}.cube.I.sd.fits"
        )

        if tp_startmodel:
            if tp_cube is None:
                print(f"[ERROR] No TP cube found for {obs_id}. Skipping.")
                continue

        twelve_m_data = glob.glob(f"{ACES_DATA / (prefix + twelve_m_mous_id) / 'calibrated/working'}/*.ms")
        seven_m_data = glob.glob(f"{ACES_DATA / (prefix + seven_m_mous_id) / 'calibrated/working'}/*.ms")

        twelve_m_visfiles = split_visfiles(twelve_m_data, obs_dir, line, array='12m')
        seven_m_visfiles = split_visfiles(seven_m_data, obs_dir, line, array='7m')

        do_mstransform(twelve_m_visfiles, obs_dir, restfreq, v_start, v_width, nchan)
        do_mstransform(seven_m_visfiles, obs_dir, restfreq, v_start, v_width, nchan)

        twelve_m_visfiles_mstransform = glob.glob(str(obs_dir) + f'/*split.SPW{LINE_SPWS[line]["mol_12m_spw"]}.mstransform')
        seven_m_visfiles_mstransform = glob.glob(str(obs_dir) + f'/*split.SPW{LINE_SPWS[line]["mol_7m_spw"]}.mstransform')

        concatvis = f"{obs_dir}/Sgr_A_st_{obs_id}.{line}.12m7m.concat.ms"
        if not Path(concatvis).exists() and seven_m_visfiles_mstransform and twelve_m_visfiles_mstransform:
            concat(vis=seven_m_visfiles_mstransform + twelve_m_visfiles_mstransform, concatvis=concatvis)

        do_clean(obs_dir, obs_id, tp_cube, concatvis, line, deep_clean, relaxed_masking, tp_startmodel, overwrite_clean)

# TODO: Add cleanup functionality to remove intermediate files
def main():
    aces_rootdir = os.getenv('ACES_ROOTDIR')
    aces_workdir = os.getenv('ACES_WORKDIR')
    aces_data = os.getenv('ACES_DATA')
    line_table = 'aces/data/tables/linelist.csv'

    if not all([aces_rootdir, aces_workdir, aces_data]):
        raise EnvironmentError("One or more required environment variables are not set.")

    init(aces_rootdir, aces_workdir, aces_data, line_table)

    # Set up the relevant parameters
    lines = ['hnco43']
    v_start = '-85 km/s'
    v_width = '0.125 km/s'
    nchan = 120
    region = ['t']
    deep_clean = True
    relaxed_masking = True
    tp_startmodel = True
    overwrite_clean = True

    for line in tqdm(lines, desc='LINES'):
        restfreq = LINE_SPWS[process_string(line)]['restfreq']
        do_joint_deconvolution(region, line, restfreq, v_start, v_width, nchan, deep_clean, relaxed_masking, tp_startmodel, overwrite_clean)

if __name__ == "__main__":
    main()