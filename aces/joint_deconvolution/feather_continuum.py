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

# MUSTANG-2 (minkasi remap of all CMZ data, calibrated against ACT+Planck) feathered with the
# ACT+Planck 90 GHz coadd: K at 90 GHz, all angular scales; the beam (10.7" core, 26% of the
# area, + 30" shoulder 26% + 115" error beam 48%) is deconvolved to the 10.7" core
MUSTANG_ACT_FEATHERED = '/orange/adamginsburg/cmz/MUSTANG/act_feather/mustang_act_feather_3c30_K.fits'
# the TENS MUSTANG-2 map of Sgr B2 (Jy/beam, no large scales: needs combine_with_planck=True)
TENS_DIR = '/orange/adamginsburg/ACES/TENS/'
TENS_MUSTANG = f'{TENS_DIR}/SgrB2_flux_cut_dt6_filtp05to49_noShift_final_map.fits'


def feather_aces_with_mustang(output_filename=None, use_cached=True, niter=20,
                              input_filename=None, remake_mustang_planck=True,
                              feather_beam=None, lowresscalefactor=1.0,
                              mustang_filename=MUSTANG_ACT_FEATHERED,
                              combine_with_planck=None):
    """
    Feather ACES continuum data with MUSTANG data.

    The default single-dish map is ``MUSTANG_ACT_FEATHERED``: the MUSTANG-2 CMZ map, calibrated
    against ACT+Planck (Guan et al. 2021) and feathered with it, in K at 90 GHz.  ACES (93.7 GHz) is
    converted to K at 97.3 GHz, so for equal flux scales the single-dish / ACES ratio in K is 1.02 for
    dust (alpha = 3.5) and 1.17 for free-free.

    Parameters
    ----------
    output_filename : str, optional
        The filename to save the feathered data to. If None, will use default path.
    use_cached : bool, optional
        If True, use cached feathered data if it exists.
    niter : int, optional
        Number of binary erosion iterations applied to the ACES footprint edge.
    input_filename : str, optional
        The ACES 12m continuum mosaic to feather.  If None, use
        ``12m_continuum_commonbeam_circular_reimaged_mosaic.fits``.
    remake_mustang_planck : bool, optional
        If False and the MUSTANG+Planck combined map already exists, use it
        instead of regenerating it (only with ``combine_with_planck``).
    feather_beam : `~astropy.units.Quantity`, optional
        Smooth the aligned single-dish map to this circular FWHM before
        feathering, so that the interferometer data supply all scales smaller
        than it.  None keeps the 15" beam of the Planck-combined map and
        smooths any other map to 15": ``feather_simple`` passes the
        single-dish map at all spatial frequencies, and the pixel noise of an
        unsmoothed minkasi map (white down to the 1" pixels) raises the noise
        of the ACES mosaic on <8" scales 6-16x (1.1-1.7 mK -> 7-25 mK).  At
        Sgr B2 M the TENS map holds about half the flux of MGPS90 (and of the
        ALMA data) at 20", and feathering it at 15" digs bowls of -2.5 to
        -5.4 K around M.
    lowresscalefactor : float, optional
        Scale factor for the single-dish map, passed to ``feather_simple``.
    mustang_filename : str, optional
        The single-dish map.  ``TENS_MUSTANG`` reproduces the earlier feathered
        mosaics (with ``combine_with_planck=True``).
    combine_with_planck : bool, optional
        Smooth the (Jy/beam, 10" beam) map to 15" and feather it with Planck
        100 GHz before use.  None: True only for ``TENS_MUSTANG``.

    Returns
    -------
    rslt : fits.PrimaryHDU
        The feathered data.
    """
    default_input_filename = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits'
    if input_filename is None:
        input_filename = default_input_filename
    if output_filename is None:
        output_filename = input_filename.replace('.fits', '_MUSTANGfeathered.fits')
        print(f"Set output_filename to {output_filename}")

    if use_cached and os.path.exists(output_filename):
        return fits.open(output_filename)[0]

    if combine_with_planck is None:
        combine_with_planck = os.path.realpath(mustang_filename) == os.path.realpath(TENS_MUSTANG)
    mustang_reffreq = 87.85e9

    # Load ACES continuum image
    ACES = fits.open(input_filename)
    ACESheader = ACES[0].header
    ACESwcs = WCS(ACES[0].header)
    ACESheader['BUNIT'] = 'K'
    ACESbeam = Beam.from_fits_header(ACESheader)
    ACESdata = ACES[0].data

    if combine_with_planck:
        # Process MUSTANG data with Planck
        mustang_beam_fwhm = 10 * u.arcsec
        mustang_central_frequency = 91.5 * u.GHz
        mustangbeam = Beam(mustang_beam_fwhm)

        # MUSTANG data and the MUSTANG+Planck product
        fn = mustang_filename
        mustang_plus_planck_fn = mustang_filename.replace('.fits', '_PlanckCombined.fits')

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

        if remake_mustang_planck or not os.path.exists(mustang_plus_planck_fn):
            # Smooth MUSTANG data to target beam (MUSTANG+Planck is regenerated unless cached and remake_mustang_planck=False)
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
    else:
        mustang_plus_planck_fn = mustang_filename

    # Load the single-dish map for feathering with ACES
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
    aces_smooth_fn = input_filename.replace('.fits', '_smoothed_to_mustang.fits')

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
        # MUSTANG appears displaced by (+xoff, +yoff) pixels, so move it by (-xoff, -yoff):
        # CRPIX += offset translates the image exactly; changing CRVAL2 of a CAR map with
        # CRVAL2 = 0 makes an oblique projection (invalid for LATPOLE = 90)
        mustang_planck_image[0].header['SHIFTED'] = (-xoff, -yoff)
        mustang_planck_image[0].header['CRPIX1'] += xoff
        mustang_planck_image[0].header['CRPIX2'] += yoff

    if feather_beam is None and not combine_with_planck:
        feather_beam = 15 * u.arcsec
    if feather_beam is not None:
        sdbeam = Beam.from_fits_header(mustang_planck_image[0].header)
        newbeam = Beam(feather_beam)
        mustang_planck_image[0].data = convolve_fft(mustang_planck_image[0].data,
                                                    newbeam.deconvolve(sdbeam).as_kernel(pixscale_mustang),
                                                    allow_huge=True)
        mustang_planck_image[0].header.update(newbeam.to_header_keywords())
        print(f"Smoothed the single-dish map from {sdbeam.major.to(u.arcsec):0.1f} to {newbeam.major.to(u.arcsec):0.1f} before feathering")

    # Reload ACES to make sure we don't multiply by jtok_aces twice
    ACES = fits.open(input_filename)
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
        lowresscalefactor=lowresscalefactor,
        return_hdu=True
    )

    # Save the feathered result
    rslt.writeto(output_filename, overwrite=True)
    print(f"Saved to {output_filename}")

    # Create a link to TENS directory if needed (only for the default product)
    if input_filename == default_input_filename:
        try:
            tens_output = f'{TENS_DIR}/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'
            os.link(output_filename, tens_output)
        except FileExistsError:
            pass  # all is good

    return rslt


if __name__ == "__main__":
    print("Running feather_aces_with_mustang...")
    rslt = feather_aces_with_mustang(use_cached=False, niter=30)
    print(rslt)