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


def get_largest_bmaj_bmin(files):
    """
    Get the largest BMAJ and BMIN values from a list of FITS files.
    """
    hdu_list = [fits.open(file)[0] for file in files]

    # Initialize largest values
    largest_bmaj, largest_bmin = 0, 0

    # Loop over HDUs
    for hdu in hdu_list:
        header = hdu.header

        # Check if BMAJ and BMIN exist in the header
        if 'BMAJ' in header and 'BMIN' in header:
            bmaj = header['BMAJ']
            bmin = header['BMIN']

            # Update largest values
            if bmaj > largest_bmaj:
                largest_bmaj = bmaj

            if bmin > largest_bmin:
                largest_bmin = bmin

    return largest_bmaj, largest_bmin


def smooth_hdus_spectral_cube(files, largest_bmaj, largest_bmin):
    """
    Smooth a set of cubes to the largest beam size.
    Note that we are only using BMAJ here so that the beam is circular.
    """
    target_beam = Beam(major=largest_bmaj*u.deg, minor=largest_bmaj*u.deg, pa=0*u.deg)
    hdu_list = [fits.open(file)[0] for file in files]
    for i in range(len(files)):

        hdu = hdu_list[i]
        cube = SpectralCube.read(hdu)
        cube.allow_huge_operations=True
        cube_beam_deg = cube.beam.major.to(u.deg).value
        if cube_beam_deg < ((largest_bmaj*u.deg).value):
            cube = cube.convolve_to(target_beam)

        print('INFO Smoothed to largest beam.')

        outfile = files[i].replace('.fits', '.smoothed.fits')
        print('INFO Save smoothed files: %s '% outfile)
        cube.write(outfile, overwrite=True)

        del cube
        gc.collect()

    return()


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


def weighted_reproject_and_coadd(cube_files, weight_files, dir_tmp='./tmp/'):
    """
    Weighted reprojection and co-addition of a series of cubes and associated weight files.

    Parameters
    ----------
    cube_files : list
        List of paths to the cubes.
    weight_files : list
        List of paths to the weights cubes.
    dir_tmp : str, optional
        Directory for storing temporary files. Defaults to './tmp/'.
    """
    print("[INFO] Reprojecting and co-adding cubes and weights.")
    assert len(cube_files) == len(weight_files), "Mismatched number of cubes and weights."

    if not os.path.isdir(dir_tmp):
        os.mkdir(dir_tmp)

    primary_hdus = [fits.open(cube_file)[0] for cube_file in cube_files]
    weight_hdus = [fits.open(weight_file)[0] for weight_file in weight_files]

    n_hdus = len(primary_hdus)
    data_reproject = []

    fake_hdus = create_fake_hdus(primary_hdus, 0)
    wcs_out, shape_out = find_optimal_celestial_wcs(fake_hdus)
    hdu_out = wcs_out.to_fits()[0]
    hdu_out.data = np.ones(shape_out)
    hdu_out.writeto('%s/hdu_out.fits' %dir_tmp, overwrite=True)

    reprojected_data, reprojected_weights = [], []

    p_bar = tqdm(range(n_hdus*2))
    p_bar.refresh()
    for i in range(n_hdus):

        print("[INFO] Processing primary_hdu[%i]" %i)

        cube = SpectralCube.read(primary_hdus[i])
        cube.allow_huge_operations=True
        cube.write('%s/cube.fits' %dir_tmp, overwrite=True)

        if i == 0:
            keys = ['CUNIT3', 'CTYPE3', 'CRPIX3', 'CDELT3',
            'CRVAL3', 'WCSAXES', 'SPECSYS', 'RESTFRQ',
            'BUNIT', 'BMAJ', 'BMIN', 'BPA']
            for key in keys:
                hdu_out.header[key] = cube.header[key]

        regrid_fits_to_template('%s/cube.fits' %dir_tmp,
                                '%s/hdu_out.fits' %dir_tmp,
                                '%s/cube_regrid_%i.fits' %(dir_tmp, i))

        cube_regrid = SpectralCube.read('%s/cube_regrid_%i.fits' %(dir_tmp, i))
        reprojected_data.append(cube_regrid.hdu.data)

        del cube
        del cube_regrid
        gc.collect()

        p_bar.update(1)
        p_bar.refresh()

        print("[INFO] Processing weight_hdus[%i]" %i)

        cube_weight = SpectralCube.read(weight_hdus[i])
        cube_weight.allow_huge_operations=True
        cube_weight.write('%s/cube_weight.fits' %dir_tmp, overwrite=True)

        regrid_fits_to_template('%s/cube_weight.fits' %dir_tmp,
                                '%s/hdu_out.fits' %dir_tmp,
                                '%s/cube_weight_regrid_%i.fits' %(dir_tmp, i))

        cube_weight_regrid = SpectralCube.read('%s/cube_weight_regrid_%i.fits' %(dir_tmp, i))
        reprojected_weights.append(cube_weight_regrid.hdu.data)

        del cube_weight
        del cube_weight_regrid
        gc.collect()
        p_bar.update(1)
        p_bar.refresh()

    weighted_data = np.array(reprojected_data) * np.array(reprojected_weights)
    data_reproject = np.nansum(weighted_data, axis=0) / np.nansum(reprojected_weights, axis=0)

    hdu_reproject = fits.PrimaryHDU(data_reproject, hdu_out.header)

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

    largest_bmaj, largest_bmin = get_largest_bmaj_bmin(cube_files)
    smooth_hdus_spectral_cube(cube_files, largest_bmaj, largest_bmin)
    cube_files = [filename.replace('.fits', '.smoothed.fits') for filename in cube_files]

    print("[INFO] Creating weighted mosaic for TP+7m+12m cubes.")
    mosaic_hdu = weighted_reproject_and_coadd(cube_files, weight_files)
    mosaic_hdu.writeto(outputfile, overwrite=True)
    print(f"[INFO] Created and saved weighted mosaic to {outputfile}")

    return()