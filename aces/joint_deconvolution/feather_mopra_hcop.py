import os
import glob
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.table import Table
from astropy.io import fits
from spectral_cube import SpectralCube
from casatasks import importfits, imregrid, exportfits, imhead, imsmooth, feather, imsubimage


def get_file(filename):
    files = glob.glob(filename)
    return str(files[0]) if files else None


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


def get_rest_frequency(image):
    if image.endswith('.fits'):
        with fits.open(image) as hdul:
            return hdul[0].header.get('RESTFRQ', hdul[0].header.get('RESTFREQ'))
    else:
        return imhead(image, mode='get', hdkey='restfreq')


def drop_stokes_axis(imagename):
    output = imagename + '.nostokes'
    if not os.path.exists(output):
        imsubimage(imagename=imagename, outfile=output, dropdeg=True)
    return output


def import_sd_fits(fitsimage, imagename):
    """
    Import single dish FITS files and apply main beam efficiency correction
    to convert from antenna temperature to main beam temperature.
    """
    if os.path.exists(imagename):
        shutil.rmtree(imagename)

    with fits.open(fitsimage) as hdul:
        header = hdul[0].header
        data = hdul[0].data

        bunit = header.get('BUNIT', '').upper()
        if 'K' in bunit:
            # Apply main beam efficiency correction
            # T_mb = T_a / eta_mb
            eta_mb = 0.49  # main beam efficiency
            data = data / eta_mb
            header['COMMENT'] = 'Data converted from Ta to Tmb using eta_mb = 0.49'

            temp_fits = fitsimage + '.temp.fits'
            fits.writeto(temp_fits, data, header, overwrite=True)
            fitsimage = temp_fits

    importfits(fitsimage=fitsimage, imagename=imagename, overwrite=True)

    if fitsimage.endswith('.temp.fits'):
        os.remove(fitsimage)

    return imagename


def import_interferometer_fits(fitsimage, imagename):
    if not os.path.exists(imagename):
        importfits(fitsimage=fitsimage, imagename=imagename, overwrite=False)
    return imagename


def match_spectral_coverage_and_resolution(sd_file, interf_file):
    from astropy.convolution import Gaussian1DKernel, convolve

    sd_cube = SpectralCube.read(sd_file)
    interf_cube = SpectralCube.read(interf_file)

    if sd_cube.spectral_axis.unit != interf_cube.spectral_axis.unit:
        sd_cube = sd_cube.with_spectral_unit(interf_cube.spectral_axis.unit)

    spec_min = interf_cube.spectral_axis.min()
    spec_max = interf_cube.spectral_axis.max()

    sd_cube_matched = sd_cube.spectral_slab(spec_min, spec_max)

    sd_chan_width = np.abs(np.median(np.diff(sd_cube_matched.spectral_axis)))
    interf_chan_width = np.abs(np.median(np.diff(interf_cube.spectral_axis)))

    sigma_sd = sd_chan_width / np.sqrt(8 * np.log(2))
    sigma_interf = interf_chan_width / np.sqrt(8 * np.log(2))

    sigma_conv = np.sqrt(sigma_sd**2 - sigma_interf**2)

    if sigma_conv > 0:
        sigma_conv_pix = (sigma_conv / interf_chan_width).decompose().value
        kernel = Gaussian1DKernel(stddev=sigma_conv_pix)
        interf_cube_smoothed = interf_cube.spectral_smooth(kernel)
    else:
        interf_cube_smoothed = interf_cube

    new_spectral_axis = sd_cube_matched.spectral_axis

    interf_cube_regridded = interf_cube_smoothed.spectral_interpolate(new_spectral_axis)

    sd_output = sd_file.replace('.fits', '.matched.fits')
    interf_output = interf_file.replace('.fits', '.regridded.fits')

    sd_cube_matched.write(sd_output, overwrite=True)
    interf_cube_regridded.write(interf_output, overwrite=True)

    return sd_output, interf_output


def prepare_sd_cube(sd_cube, base_sd_name):
    sd_cube_im = base_sd_name + '.image'
    sd_cube_eq = base_sd_name + '.j2000'
    sd_cube_eq_fits = base_sd_name + '.j2000.fits'

    if not os.path.exists(sd_cube_im):
        import_sd_fits(sd_cube, sd_cube_im)

    if not os.path.exists(sd_cube_eq):
        imregrid(imagename=sd_cube_im, template='J2000', output=sd_cube_eq)

    if not os.path.exists(sd_cube_eq_fits):
        exportfits(imagename=sd_cube_eq, fitsimage=sd_cube_eq_fits,
                   dropdeg=True, overwrite=True)
 
    return sd_cube_eq_fits


def prepare_and_feather(obs_dir, obs_id, sd_cube_eq_fits, seven_m_cube, twelve_m_cube, MOLECULE, cleanup=True):
    """Prepare and feather data in two steps: First SD+7M, then (SD+7M)+12M."""
    OUTPUT_SD_7M = obs_dir / f'Sgr_A_st_{obs_id}.SD_7M_feather.{MOLECULE}.image.statcont.contsub'
    OUTPUT_ALL = obs_dir / f'Sgr_A_st_{obs_id}.SD_7M_12M_feather.{MOLECULE}.image.statcont.contsub'

    if os.path.isfile(str(OUTPUT_ALL) + '.fits'):
        print(f"Final output for region {obs_id} already exists, skipping...")
        return True

    print(f"Processing region {obs_id}...")

    # Step 1: SD+7M feathering
    seven_m_head = imhead(seven_m_cube)
    if 'perplanebeams' in seven_m_head:
        commonbeam = seven_m_cube.replace('.fits', '.image.commonbeam')
        if not os.path.isdir(commonbeam):
            imsmooth(imagename=seven_m_cube, kernel='commonbeam',
                     targetres=True, outfile=commonbeam)
        if not os.path.isfile(commonbeam + '.fits'):
            exportfits(imagename=commonbeam, fitsimage=commonbeam + '.fits', dropdeg=True)
        seven_m_cube = commonbeam + '.fits'

    sd_cube_matched, seven_m_cube_regridded = match_spectral_coverage_and_resolution(
        sd_cube_eq_fits, seven_m_cube
    )

    sd_cube_matched_im = import_sd_fits(
        sd_cube_matched,
        sd_cube_matched.replace('.fits', '.image')
    )

    seven_m_cube_im = import_interferometer_fits(
        seven_m_cube_regridded,
        seven_m_cube_regridded.replace('.fits', '.image')
    )

    seven_m_cube_im = drop_stokes_axis(seven_m_cube_im)

    if not os.path.exists(str(OUTPUT_SD_7M)):
        feather(imagename=str(OUTPUT_SD_7M), highres=seven_m_cube_im,
                lowres=sd_cube_matched_im, sdfactor=1.0)

    if not os.path.isfile(str(OUTPUT_SD_7M) + '.fits'):
        exportfits(imagename=str(OUTPUT_SD_7M), fitsimage=str(OUTPUT_SD_7M) + '.fits',
                   dropdeg=True)

    # Clean up temporary files from first feathering step
    if cleanup:
        for temp_file in [sd_cube_matched, seven_m_cube_regridded]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        if os.path.exists(sd_cube_matched_im):
            shutil.rmtree(sd_cube_matched_im)

    # Step 2: (SD+7M)+12M feathering
    if twelve_m_cube:
        twelve_m_head = imhead(twelve_m_cube)
        if 'perplanebeams' in twelve_m_head:
            commonbeam = twelve_m_cube.replace('.fits', '.image.commonbeam')
            if not os.path.isdir(commonbeam):
                imsmooth(imagename=twelve_m_cube, kernel='commonbeam',
                         targetres=True, outfile=commonbeam)
            if not os.path.isfile(commonbeam + '.fits'):
                exportfits(imagename=commonbeam, fitsimage=commonbeam + '.fits', dropdeg=True)
            twelve_m_cube = commonbeam + '.fits'

        # Prepare 12m data for feathering with SD+7M data
        _, twelve_m_regridded = match_spectral_coverage_and_resolution(
            str(OUTPUT_SD_7M) + '.fits', twelve_m_cube
        )

        twelve_m_cube_im = import_interferometer_fits(
            twelve_m_regridded,
            twelve_m_regridded.replace('.fits', '.image')
        )

        twelve_m_cube_im = drop_stokes_axis(twelve_m_cube_im)

        if not os.path.exists(str(OUTPUT_ALL)):
            feather(imagename=str(OUTPUT_ALL), highres=twelve_m_cube_im,
                    lowres=str(OUTPUT_SD_7M), sdfactor=1.0)

        if not os.path.isfile(str(OUTPUT_ALL) + '.fits'):
            exportfits(imagename=str(OUTPUT_ALL), fitsimage=str(OUTPUT_ALL) + '.fits',
                       dropdeg=True)

        # Clean up temporary files from second feathering step
        if cleanup:
            for temp_file in [twelve_m_regridded]:
                if os.path.exists(temp_file):
                    os.remove(temp_file)

    return True


def create_feathercubes(ACES_WORKDIR, ACES_DATA, ACES_ROOTDIR, line_spws, MOLECULE):
    sb_names = pd.read_csv(ACES_ROOTDIR / 'aces/data/tables/aces_SB_uids.csv')
    sd_cube = "CMZ_3mm_HCO.fits"
    generic_name = '.Sgr_A_star_sci.spw'
    prefix = 'member.uid___A001_'

    base_sd_name = Path(sd_cube).stem
    sd_cube_eq = prepare_sd_cube(sd_cube, base_sd_name)

    for i, row in sb_names.iterrows():
        obs_id = row['Obs ID']
        obs_dir = ACES_WORKDIR / f'Sgr_A_st_{obs_id}'
        obs_dir.mkdir(exist_ok=True)

        seven_m_cube = get_file(
            f"{ACES_DATA / (prefix + row['7m MOUS ID']) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_7m_spw']}.cube.I.iter1.image.pbcor.statcont.contsub.fits"
        )

        twelve_m_cube = get_file(
            f"{ACES_DATA / (prefix + row['12m MOUS ID']) / 'calibrated/working'}/*{generic_name}{line_spws[MOLECULE]['mol_12m_spw']}.cube.I.iter1.image.pbcor.statcont.contsub.fits"
        )

        if seven_m_cube and twelve_m_cube:
            prepare_and_feather(obs_dir, obs_id, sd_cube_eq, seven_m_cube, twelve_m_cube, MOLECULE)


if __name__ == "__main__":
    ACES_ROOTDIR = Path(os.getenv('ACES_ROOTDIR'))
    ACES_WORKDIR = Path(os.getenv('ACES_DATA'))
    ACES_DATA = Path(os.getenv('ACES_DATA'))
    LINE_TABLE = (ACES_ROOTDIR / 'aces/data/tables/linelist.csv')
    LINES = ['hco+10']

    for LINE in LINES:
        LINE_SPWS = read_table(LINE_TABLE)
        create_feathercubes(ACES_WORKDIR, ACES_DATA, ACES_ROOTDIR, LINE_SPWS, LINE)