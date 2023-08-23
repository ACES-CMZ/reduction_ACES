import os
import gc
import warnings
import numpy as np
from tqdm import tqdm
from radio_beam import Beam
from astropy.io import fits
from astropy import units as u
from spectral_cube import SpectralCube
from reproject.mosaicking import find_optimal_celestial_wcs
from casatasks import exportfits, importfits, imregrid, rmtables
warnings.filterwarnings('ignore')


def create_fake_hdus(hdus, j):
    """
    Create a set of fake HDUs.
    This is a bit of a workaround so that we can use reproject/find_optimal_celestial_wcs with 3D data.
    """
    print("[INFO] Creating fake HDUs for reprojecting.")
    fake_hdus = []
    for i in range(len(hdus)):
        data = hdus[i].data.copy()
        header = hdus[i].header.copy()
        data = data[j]
        del header['*3*']
        del header['*4*']
        header['WCSAXES'] = 2
        fake_hdus.append(fits.PrimaryHDU(data, header))
    return fake_hdus


def get_common_beam(files):
    """
    Get the common beam size of a set of cubes.
    """
    cube = SpectralCube.read(files[0])
    common_beam = cube.beam

    for fn in files[1:]:
        cube = SpectralCube.read(fn)
        common_beam = common_beam.commonbeam_with(cube.beam)
    return common_beam


def smooth_to_common_beam(files, common_beam):
    """
    Smooth a set of cubes to common beam size.
    """
    hdu_list = [fits.open(file)[0] for file in files]
    for i in range(len(files)):
        if not os.path.exists(files[i].replace('.fits', '.smoothed.fits')):
            hdu = hdu_list[i]
            cube = SpectralCube.read(hdu)
            cube.allow_huge_operations = True
            cube = cube.convolve_to(common_beam)

            print('INFO Smoothed to common beam.')
            outfile = files[i].replace('.fits', '.smoothed.fits')
            print('INFO Save smoothed files: %s ' % outfile)
            cube.write(outfile, overwrite=True)


def regrid_fits_to_template(input_fits, template_fits, output_fits, overwrite=True):
    """
    Regrid a FITS file using a template image and save the regridded cube to a new FITS file.
    """
    if os.path.exists(output_fits) and not overwrite:
        print("Output file already exists. Use `overwrite=True` to overwrite it.")
        return

    # Define the names of the intermediate images
    input_image = input_fits.replace('.fits', '.tmp.image')
    template_image = template_fits.replace('.fits', '.tmp.image')
    regrid_image = input_image.replace('.image', '.tmp.regrid.image')

    # Remove any pre-existing intermediate images
    rmtables([input_image, template_image, regrid_image])

    importfits(fitsimage=input_fits, imagename=input_image, overwrite=overwrite)
    importfits(fitsimage=template_fits, imagename=template_image, overwrite=overwrite)
    imregrid(imagename=input_image, template=template_image, output=regrid_image)
    exportfits(imagename=regrid_image, fitsimage=output_fits, overwrite=overwrite)

    # Clean up by removing the intermediate images
    rmtables([input_image, template_image, regrid_image])


def weighted_reproject_and_coadd(cube_files, weight_files, dir_tmp='./tmp/', overwrite_dir_tmp=False):
    '''
    Function to reproject and coadd the cubes and weights.

    Inputs:
    - cube_files: a list of paths to the cube files
    - weight_files: a list of paths to the weight files

    Outputs:
    - a HDU object representing the reprojected and coadded data
    '''
    tqdm.write("[INFO] Reprojecting and co-adding cubes and weights.")
    assert len(cube_files) == len(weight_files), "Mismatched number of cubes and weights."

    if overwrite_dir_tmp:
        os.system('rm -rf %s' % dir_tmp)

    if not os.path.isdir(dir_tmp):
        os.mkdir(dir_tmp)

    if overwrite_dir_tmp:

        tqdm.write("Processing fake hdu data")

        primary_hdus = [fits.open(cube_file)[0] for cube_file in cube_files]
        weight_hdus = [fits.open(weight_file)[0] for weight_file in weight_files]

        fake_hdus = create_fake_hdus(primary_hdus, 0)
        wcs_out, shape_out = find_optimal_celestial_wcs(fake_hdus)
        hdu_out = wcs_out.to_fits()[0]
        hdu_out.data = np.ones(shape_out)
        hdu_out.writeto('%s/hdu_out.fits' % dir_tmp, overwrite=True)
    else:
        hdu_out = fits.open('%s/hdu_out.fits' % dir_tmp)[0]

        cube = fits.open('%s/cube.fits' % dir_tmp)[0]
        keys = ['CUNIT3', 'CTYPE3', 'CRPIX3', 'CDELT3', 'CRVAL3', 'SPECSYS', 'RESTFRQ', 'BUNIT', 'BMAJ', 'BMIN', 'BPA']
        for key in keys:
            hdu_out.header[key] = cube.header[key]

        tqdm.write("Skipping processing of fake hdu data")

    n_hdus = len(cube_files)
    data_reproject = []
    reprojected_data, reprojected_weights = [], []

    p_bar = tqdm(range(n_hdus * 2))
    p_bar.refresh()
    for i in range(n_hdus):
        if os.path.isfile('%s/cube_regrid_%i.fits' % (dir_tmp, i)):
            tqdm.write("[INFO] Exists, not processing primary_hdu[%i]" % i)
            cube_regrid = fits.open('%s/cube_regrid_%i.fits' % (dir_tmp, i))[0]
        else:
            tqdm.write("[INFO] Processing primary_hdu[%i]" % i)
            cube = primary_hdus[i]
            cube.writeto('%s/cube.fits' % dir_tmp, overwrite=True)

            if i == 0:
                keys = ['CUNIT3', 'CTYPE3', 'CRPIX3', 'CDELT3', 'CRVAL3', 'SPECSYS', 'RESTFRQ', 'BUNIT', 'BMAJ', 'BMIN', 'BPA']
                for key in keys:
                    hdu_out.header[key] = cube.header[key]

            regrid_fits_to_template('%s/cube.fits' % dir_tmp,
                                    '%s/hdu_out.fits' % dir_tmp,
                                    '%s/cube_regrid_%i.fits' % (dir_tmp, i))

            cube_regrid = fits.open('%s/cube_regrid_%i.fits' % (dir_tmp, i))[0]
            del cube

        reprojected_data.append(cube_regrid.data)
        del cube_regrid
        gc.collect()

        p_bar.update(1)
        p_bar.refresh()

        if os.path.isfile('%s/cube_weight_regrid_%i.fits' % (dir_tmp, i)):

            tqdm.write("Exists, not processing weight_hdus[%i]" % i)
            cube_weight_regrid = fits.open('%s/cube_weight_regrid_%i.fits' % (dir_tmp, i))[0]
        else:
            tqdm.write("[INFO] Processing weight_hdus[%i]" % i)
            cube_weight = weight_hdus[i]
            cube_weight.writeto('%s/cube_weight.fits' % dir_tmp, overwrite=True)

            regrid_fits_to_template('%s/cube_weight.fits' % dir_tmp,
                                    '%s/hdu_out.fits' % dir_tmp,
                                    '%s/cube_weight_regrid_%i.fits' % (dir_tmp, i))

            cube_weight_regrid = fits.open('%s/cube_weight_regrid_%i.fits' % (dir_tmp, i))[0]
            del cube_weight

        reprojected_weights.append(cube_weight_regrid.data)
        del cube_weight_regrid
        gc.collect()
        p_bar.update(1)
        p_bar.refresh()

    tqdm.write('[INFO] Creating weighted_data')
    weighted_data = np.array(reprojected_data, dtype=np.float32) * np.array(reprojected_weights, dtype=np.float32)
    del reprojected_data
    gc.collect()

    tqdm.write('[INFO] Summing weighted_data --> weighted_data_sum')
    weighted_data_sum = np.nansum(weighted_data, axis=0)
    del weighted_data
    gc.collect()

    tqdm.write('[INFO] Summing reprojected_weights --> reprojected_weights_sum')
    reprojected_weights_sum = np.nansum(reprojected_weights, axis=0)
    del reprojected_weights
    gc.collect()

    tqdm.write('[INFO] Creating data_reproject')
    data_reproject = weighted_data_sum / reprojected_weights_sum
    del weighted_data_sum
    del reprojected_weights_sum
    gc.collect()

    tqdm.write('[INFO] Creating hdu_reproject')
    hdu_reproject = fits.PrimaryHDU(data_reproject, hdu_out.header)
    del data_reproject
    gc.collect()

    return hdu_reproject


def create_weighted_mosaic(ACES_WORKDIR, START_VELOCITY, END_VELOCITY, VEL_RES, MOLECULE, process_12M=True):
    """
    Create a weighted mosaic of input cubes and save it to a new FITS file.

    The function first smoothes all the cubes to the same beam size, then performs a weighted reprojection and
    co-addition of the cubes using the associated weight files.

    Parameters
    ----------
    ACES_WORKDIR : Path
        Path object representing the working directory.
    START_VELOCITY : float
        The start of the velocity range in km/s.
    END_VELOCITY : float
        The end of the velocity range in km/s.
    VEL_RES : float
        The velocity resolution in km/s.
    MOLECULE : str
        The molecule being processed (e.g. 'HNCO').
    process_12M : bool, optional
        If True, the filenames are specified to include 12m data. Default is True.
    """
    if process_12M:
        cube_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.TP_7M_12M_feather_all.{MOLECULE}.image.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.fits')]
        weight_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.12M.{MOLECULE}.image.weight.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.fits')]
        outputfile = ACES_WORKDIR / f'{MOLECULE}.TP_7M_12M_weighted_mosaic.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.fits'
    else:
        cube_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.TP_7M_feather_all.{MOLECULE}.image.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.fits')]
        weight_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.7M.{MOLECULE}.image.weight.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.fits')]
        outputfile = ACES_WORKDIR / f'{MOLECULE}.TP_7M_weighted_mosaic.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.fits'

    common_beam = get_common_beam(cube_files)
    smooth_to_common_beam(cube_files, common_beam)
    cube_files = [filename.replace('.fits', '.smoothed.fits') for filename in cube_files]

    print("[INFO] Creating weighted mosaic for TP+7m+12m cubes.")
    mosaic_hdu = weighted_reproject_and_coadd(cube_files, weight_files)
    mosaic_hdu.writeto(outputfile, overwrite=True)
    print(f"[INFO] Created and saved weighted mosaic to {outputfile}")