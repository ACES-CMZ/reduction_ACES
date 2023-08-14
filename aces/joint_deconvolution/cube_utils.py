import warnings
import numpy as np
from pathlib import Path
import astropy.units as u
from astropy.io import fits
from spectral_cube import SpectralCube
from casatasks import exportfits, importfits, rmtables, imrebin
warnings.filterwarnings('ignore')


def convert_to_float32(hdu):
    """
    Does this one really need an explanation? Or even a separate function for that matter ...
    """
    hdu.data = hdu.data.astype('float32')
    return hdu


def crop_cube_velocity_range(fits_file, rest_FREQUENCY, v_start, v_end, v_res=None):
    """
    Crop cube to a specified velocity range and optionally regrid the cropped cube
    to a specified velocity resolution.

    If the cropped cube only contains one spectral channel or less, the function returns None. The function
    also returns the cube in the minimal subcube that encompasses all the non-blank pixels.

    Parameters
    ----------
    fits_file : str
        Path to the FITS file containing the cube.
    rest_FREQUENCY : float
        Rest frequency value in GHz to convert the cube's spectral axis to velocity.
    v_start : float
        The start of the velocity range in km/s.
    v_end : float
        The end of the velocity range in km/s.
    v_res : float, optional
        The velocity resolution in km/s for regridding. If not provided, no regridding is performed.
    """
    cube = SpectralCube.read(fits_file)
    cube.allow_huge_operations = True
    cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=rest_FREQUENCY * u.GHz)

    vrange = [v_start * u.km / u.s, v_end * u.km / u.s]
    cropped_cube = cube.spectral_slab(*vrange)

    if cropped_cube.shape[0] <= 1:
        return None

    if v_res is not None:
        cropped_cube = regrid_cube(cropped_cube, v_start, v_end, v_res)

    cropped_cube = cropped_cube.minimal_subcube()

    hdu = cropped_cube.hdu
    hdu = fits.PrimaryHDU(hdu.data, hdu.header)
    hdu = convert_to_float32(hdu)
    return hdu


def regrid_cube(cube, v_start, v_end, v_res):
    new_velocity = np.arange(v_start, v_end + v_res, v_res) * u.km / u.s
    new_cube = cube.spectral_interpolate(new_velocity, suppress_smooth_warning=True)
    return new_cube


def crop_cubes(ACES_WORKDIR, START_VELOCITY, END_VELOCITY, VEL_RES, line_spws, line, process_12M=True):
    """
    Crop and optionally regrid a series of cubes and their associated weight files.
    The function looks for cube files and weight files in the provided working directory and its subdirectories.
    The cropped and regridded cubes and weights are written to new FITS files.

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
    line_spws : dict
        Dictionary mapping molecule identifiers to spectral window information.
    line : str
        The molecule being processed (e.g. 'HNCO').
    process_12M : bool, optional
        If True, the filenames are specified to include 12m data. Default is True.
    """
    def crop_cube(input_file, output_file):
        if not Path(output_file).exists():
            hdu = crop_cube_velocity_range(input_file, line_spws[line]['restfreq'], v_start=START_VELOCITY, v_end=END_VELOCITY, v_res=VEL_RES)
            if hdu is None:
                print(f'Could not crop {input_file}')
            else:
                hdu.writeto(output_file)

    if process_12M:
        cube_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.TP_7M_12M_feather_all.{line}.image.fits')]
        weight_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.12M.{line}.image.weight.fits')]
    else:
        cube_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.TP_7M_feather_all.{line}.image.fits')]
        weight_files = [str(x) for x in ACES_WORKDIR.glob(f'**/*.7M.{line}.image.weight.fits')]

    for cube_file, weight_file in zip(sorted(cube_files), sorted(weight_files)):
        outputfile_cube = cube_file.replace('.fits', f'.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.fits')
        outputfile_weight = weight_file.replace('.fits', f'.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.fits')

        crop_cube(cube_file, outputfile_cube)
        crop_cube(weight_file, outputfile_weight)


def cubeconvert_K_kms(ACES_WORKDIR, MOLECULE, START_VELOCITY, END_VELOCITY, VEL_RES, process_12M=True):
    """
    Convert the units of a cube to Kelvin and velocity (km/s), and save the converted cube to a new FITS file.

    Parameters
    ----------
    ACES_WORKDIR : Path
        Path object representing the working directory.
    MOLECULE : str
        The molecule being processed (e.g. 'HNCO').
    START_VELOCITY : float
        The start of the velocity range in km/s.
    END_VELOCITY : float
        The end of the velocity range in km/s.
    VEL_RES : float
        The velocity resolution in km/s.
    process_12M : bool, optional
        If True, the filename is specified to include 12m data. Defaults to True.
    """
    if process_12M:
        fits_file = f'{ACES_WORKDIR}/{MOLECULE}.TP_7M_12M_weighted_mosaic.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.fits'
    else:
        fits_file = f'{ACES_WORKDIR}/{MOLECULE}.TP_7M_weighted_mosaic.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.fits'

    cube = SpectralCube.read(fits_file)
    cube.allow_huge_operations = True

    cube = cube.with_spectral_unit(u.km / u.s, velocity_convention='radio')
    cube = cube.minimal_subcube()
    cube = cube.to(u.K)
    hdu = cube.hdu
    hdu = fits.PrimaryHDU(hdu.data, hdu.header)
    hdu = convert_to_float32(hdu)

    outputfile = fits_file.replace('.fits', '.K.kms.fits')
    print(f"[INFO] Created and saved weighted mosaic to {outputfile}")
    hdu.writeto(outputfile, overwrite=True)

    return None


def rebin(ACES_WORKDIR, MOLECULE, START_VELOCITY, END_VELOCITY, VEL_RES, REBIN_FACTOR, process_12M=True, overwrite=True):
    """
    Rebin a cube by a given factor along the spatial axes and save the rebinned cube to a new FITS file.

    Parameters
    ----------
    ACES_WORKDIR : Path
        Path object representing the working directory.
    MOLECULE : str
        The molecule being processed (e.g. 'HNCO')
    START_VELOCITY : float
        The start of the velocity range in km/s.
    END_VELOCITY : float
        The end of the velocity range in km/s.
    VEL_RES : float
        The velocity resolution in km/s.
    REBIN_FACTOR : int
        The factor by which to rebin the spatial axes of the cube.
    process_12M : bool, optional
        If True, the filename is specified to include 12m data. Defaults to True.
    overwrite : bool, optional
        If True, overwrite any existing files with the same name. Defaults to True.
    """
    if process_12M:
        input_fits = f'{ACES_WORKDIR}/{MOLECULE}.TP_7M_12M_weighted_mosaic.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.K.kms.fits'
    else:
        input_fits = f'{ACES_WORKDIR}/{MOLECULE}.TP_7M_weighted_mosaic.{START_VELOCITY}_to_{END_VELOCITY}_kms.{VEL_RES}_kms_resolution.K.kms.fits'

    # Define the names of the intermediate images
    input_image = input_fits.replace('.fits', '.tmp.image')
    regrid_image = input_fits.replace('.fits', '.tmp.regrid.image')
    output_fits = input_fits.replace('.fits', '.rebin.fits')

    # Remove any pre-existing intermediate images
    rmtables([input_image, regrid_image])

    importfits(fitsimage=input_fits, imagename=input_image, overwrite=overwrite)
    imrebin(imagename=input_image, outfile=regrid_image, factor=[REBIN_FACTOR, REBIN_FACTOR, 1], overwrite=True)
    exportfits(imagename=regrid_image, fitsimage=output_fits, overwrite=overwrite, velocity=True)
    rmtables([input_image, regrid_image])