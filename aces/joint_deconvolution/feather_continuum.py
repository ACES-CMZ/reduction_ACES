"""
Script to feather ACES continuum data with MUSTANG data.
This script extracts functionality from the uvcombine_mustang notebook.
"""

import os
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy import wcs
from astropy.convolution import convolve_fft, Gaussian2DKernel
from radio_beam import Beam
from uvcombine import feather_simple
from spectral_cube import SpectralCube, Projection
from scipy import ndimage
from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord
from image_registration import chi2_shift

from aces import conf
basepath = conf.basepath

def feather_aces_with_mustang(output_filename=None, use_cached=True, niter=20):
    """
    Feather ACES continuum data with MUSTANG data.
    
    Parameters
    ----------
    output_filename : str, optional
        The filename to save the feathered data to. If None, will use default path.
    use_cached : bool, optional
        If True, use cached feathered data if it exists.
    
    Returns
    -------
    rslt : fits.PrimaryHDU
        The feathered data.
    """
    if output_filename is None:
        output_filename = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'
        print(f"Set output_filename to {output_filename}")
    
    if use_cached and os.path.exists(output_filename):
        return fits.open(output_filename)[0]
    
    # Define TENS directory
    TENS_DIR = '/orange/adamginsburg/ACES/TENS/'
    mustang_name = 'SgrB2_flux_cut_dt6_filtp05to49_noShift_final_map.fits'
    mustang_reffreq = 87.85e9
    
    # Load ACES continuum image
    ACES = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits')
    ACESheader = ACES[0].header
    ACESwcs = WCS(ACES[0].header)
    ACESheader['BUNIT'] = 'K'
    ACESbeam = Beam.from_fits_header(ACESheader)
    ACESdata = ACES[0].data

    # Process MUSTANG data with Planck
    mustang_beam_fwhm = 10 * u.arcsec
    mustang_central_frequency = 91.5 * u.GHz
    mustangbeam = Beam(mustang_beam_fwhm)
    
    # Path to MUSTANG data
    fn = f"{TENS_DIR}/{mustang_name}"
    mustang_plus_planck_fn = f"{TENS_DIR}/{mustang_name.replace('.fits','')}_PlanckCombined.fits"
    
    # Open MUSTANG data and prepare WCS
    fh = fits.open(fn)
    mustangwcs = wcs.WCS(fh[0].header)
    center = SkyCoord(fh[0].header['CRVAL1'], fh[0].header['CRVAL2'],
                      frame=wcs.utils.wcs_to_celestial_frame(mustangwcs),
                      unit=(u.deg, u.deg))
    
    # Get or create Planck data
    planckfn = f'{TENS_DIR}/PLCKI_l0.0b0.0_100GHz.fits'
    planckfn_scaled = planckfn.replace(".fits", "scaled_to_92.44GHz.fits")
    reffreq_hz = 104.225e9
    
    if os.path.exists(planckfn):
        planck_image = fits.open(planckfn)[0]
        planck_image.header['REFFREQ'] = reffreq_hz
        planck_image.header['BMAJ'] = 9.65/60
        planck_image.header['BMIN'] = 9.65/60
        planck_image.header['BUNIT'] = 'K'
    else:
        planck_image = SkyView.get_images(center, 'Planck 100 I', pixels=600)[0][0]
        planck_image.header['REFFREQ'] = reffreq_hz
        planck_image.header['BMAJ'] = 9.65/60
        planck_image.header['BMIN'] = 9.65/60
        planck_image.header['BUNIT'] = 'K'
        planck_image.writeto(planckfn, overwrite=True)
    
    # Scale Planck data to MUSTANG frequency
    if not os.path.exists(planckfn_scaled):
        planck_image.data = planck_image.data * (reffreq_hz / mustang_reffreq)**-3
        planck_image.header['REFFREQ'] = mustang_reffreq
        planck_image.writeto(planckfn_scaled, overwrite=True)
    
    # Set MUSTANG beam and frequency information
    fh[0].header['BMAJ'] = mustang_beam_fwhm.to(u.deg).value
    fh[0].header['BMIN'] = mustang_beam_fwhm.to(u.deg).value
    fh[0].header['REFFREQ'] = mustang_central_frequency.to(u.Hz).value
    fh[0].header['BUNIT'] = 'Jy/beam'
    
    # Smooth MUSTANG data to target beam
    target_beam = Beam(15*u.arcsec)
    pixscale_mustang = (mustangwcs.proj_plane_pixel_area()**0.5).to(u.arcsec)
    kernel = target_beam.deconvolve(mustangbeam).as_kernel(pixscale_mustang)
    
    sm = convolve_fft(fh[0].data.astype('float64'), kernel).real
    fh[0].header['BMAJ'] = target_beam.major.to(u.deg).value
    fh[0].header['BMIN'] = target_beam.minor.to(u.deg).value
    
    # Convert MUSTANG data from Jy/beam to K
    jtok_mustang = mustangbeam.jtok(mustang_central_frequency)
    header_copy = fh[0].header.copy()
    header_copy['BUNIT'] = 'K'
    fh_K = fits.PrimaryHDU(sm*jtok_mustang, header=header_copy)
    
    # Feather MUSTANG+Planck
    rslt = feather_simple(fh_K, planckfn_scaled)
    
    # Save MUSTANG+Planck
    fh[0].data = rslt.real
    fh[0].header['BUNIT'] = 'K'
    hdu = fits.PrimaryHDU(data=rslt.real.astype('>f8'), header=fh[0].header.copy())
    hdu.writeto(mustang_plus_planck_fn, overwrite=True, checksum=True)
    
    # Load MUSTANG+Planck data for feathering with ACES
    mustang_planck_image = fits.open(mustang_plus_planck_fn)
    mustangdata = mustang_planck_image[0].data
    mustangwcs = WCS(mustang_planck_image[0].header)
    mustang_header = mustang_planck_image[0].header
    
    # Convert to Kelvin if not already
    mustangbeam = Beam.from_fits_header(mustang_header)
    mustang_header['REFFREQ'] = mustang_reffreq
    if mustang_header['BUNIT'].strip().upper() == 'JY/BEAM':
        mustangdata *= (1*u.Jy).to(u.K, mustangbeam.jtok_equiv(mustang_header['REFFREQ']*u.Hz)).value
        mustang_header['BUNIT'] = 'K'
    
    # Check alignment between ACES and MUSTANG
    # First smooth ACES to match MUSTANG resolution
    aces_smooth_fn = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_smoothed_to_mustang.fits'
    
    # Convert ACES data to Kelvin
    jtok_aces = (1*u.Jy).to(u.K, Beam.from_fits_header(ACESheader).jtok_equiv(97.3*u.GHz)).value
    ACESdata *= jtok_aces
    
    # Smooth ACES and reproject to MUSTANG grid
    convbeam = mustangbeam.deconvolve(ACESbeam)
    pixscale = ACESwcs.proj_plane_pixel_area()**0.5
    header = ACESheader.copy()
    header.update(mustang_header)  # replace WCS with MUSTANG wcs, but keep other metadata
    header.update(mustangbeam.to_header_keywords())
    header['SMOOTHED'] = 'True'
    
    ACES_smoothed = convolve_fft(ACESdata, convbeam.as_kernel(pixscale))
    
    import reproject
    ACES_smooth_repr, footprint = reproject.reproject_interp(
        (ACES_smoothed, ACESheader),
        mustang_header
    )
    
    ACES_smoothed_hdu = fits.PrimaryHDU(data=ACES_smooth_repr, header=header)
    ACES_smoothed_hdu.writeto(aces_smooth_fn, overwrite=True)
    
    # Align MUSTANG to ACES if needed
    xoff, yoff, _, _ = chi2_shift(ACES_smoothed_hdu.data, mustangdata, verbose=True)
    
    pixscale_mustang = (mustangwcs.proj_plane_pixel_area()**0.5).to(u.arcsec)
    print(f"MUSTANG is shifted from ACES by dl={xoff*pixscale_mustang:0.3f}, db={yoff*pixscale_mustang:0.3f}")
    
    if 'SHIFTED' not in mustang_planck_image[0].header:
        mustang_planck_image[0].header['SHIFTED'] = (-xoff, -yoff)
        mustang_planck_image[0].header['CRVAL1'] -= (xoff*pixscale_mustang).to(u.deg).value
        mustang_planck_image[0].header['CRVAL2'] -= (yoff*pixscale_mustang).to(u.deg).value
    
    # Reload ACES to make sure we don't multiply by jtok_aces twice
    ACES = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits')
    ACESheader = ACES[0].header
    ACESwcs = WCS(ACES[0].header)
    ACESheader['BUNIT'] = 'K'
    ACESbeam = Beam.from_fits_header(ACESheader)
    ACESdata = ACES[0].data
    
    # Convert ACES data to Kelvin
    jtok_aces = (1*u.Jy).to(u.K, Beam.from_fits_header(ACESheader).jtok_equiv(97.3*u.GHz)).value
    ACESdata *= jtok_aces

    # Added 2025-04-06: there are some nasty wrinkly edges that we can clean up
    aces_footprint = np.isfinite(ACESdata)
    # 4 iterations gets ride of most of the negatives, 14 gets rid of the noisy edges...
    eroded = ndimage.binary_erosion(aces_footprint, iterations=niter)
    ACESdata[~eroded] = np.nan
    ACES[0].data = ACESdata # redundant, but doing it just in case...
    print(f"Eroded ACES data to remove negative values with {niter} iterations")
    
    # Feather ACES and MUSTANG
    rslt = feather_simple(
        hires=ACES[0],
        lores=mustang_planck_image[0],
        return_hdu=True
    )
    
    # Save the feathered result
    rslt.writeto(output_filename, overwrite=True)
    print(f"Saved to {output_filename}")
    
    # Create a link to TENS directory if needed
    try:
        tens_output = f'{TENS_DIR}/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'
        os.link(output_filename, tens_output)
    except FileExistsError:
        pass  # all is good
    
    return rslt

if __name__ == "__main__":
    print("Running feather_aces_with_mustang...")
    rslt = feather_aces_with_mustang(use_cached=False)
    print(rslt)