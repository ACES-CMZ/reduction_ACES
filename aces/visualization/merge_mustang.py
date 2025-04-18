# flake8: noqa
"""
This is probably mostly done in a notebook...

it doesn't work
"""
from astropy import visualization
from astropy import wcs
import uvcombine
from astropy import coordinates
from astropy import units as u
from radio_beam import Beam
from astropy.io import fits
import numpy as np
from astropy.convolution import convolve_fft, Gaussian2DKernel

from astroquery.skyview import SkyView
from uvcombine import feather_simple, feather_plot
from astropy import wcs

from astropy import units as u
from astropy import coordinates
from astropy import wcs

from radio_beam import Beam

# based on alpha=0
mustang_reffreq = 87.85e9


def get_mustang_data_old():
    mustangfn = '/orange/adamginsburg/mgps/mgps/SgrB2/SgrB2_5pass_1_.0.2_10mJy_10mJy_w_session5_final_smooth4_PlanckCombined.fits'
    mustang_planck_image = fits.open(mustangfn)
    mustangdata = mustang_planck_image[0].data
    mustangwcs = wcs.WCS(mustang_planck_image[0].header)
    mustangheader = mustang_planck_image[0].header
    # convert to Kelvin
    mustangbeam = Beam.from_fits_header(mustangheader)

    # set reference frequency to 87.85 GHz for alpha=0
    mustangheader['REFFREQ'] = 87.85e9
    mustangdata *= (1*u.Jy).to(u.K, mustangbeam.jtok_equiv(mustangheader['REFFREQ']*u.Hz)).value
    mustangheader['BUNIT'] = 'K'

    return mustangdata, mustangheader, mustangwcs, mustangbeam


def get_ACES_data(filename='12m_continuum_commonbeam_circular_reimaged_mosaic.fits',
                  path='/orange/adamginsburg/ACES/mosaics/continuum/'):
    acesfn = f'{path}/{filename}'
    ACES = fits.open(acesfn)
    ACESheader = ACES[0].header
    ACESwcs = wcs.WCS(ACES[0].header)
    ACESheader['BUNIT'] = 'K'
    ACESbeam = Beam.from_fits_header(ACESheader)
    ACESdata = ACES[0].data

    return ACESdata, ACESheader, ACESwcs, ACESbeam


def get_mustang_data():

    mustang_beam_fwhm = 10*u.arcsec

    mgps_beam = Beam(mustang_beam_fwhm)

    #fn = 'AGBT23A_268_niter_2_25_gal.fits'
    fn = 'SgrB2_flux_cut_dt6_filtp05to49_noShift_final_map.fits'
    mustang_plus_planck_fn = f"{TENS_DIR}/{fn.replace('.fits','')}_PlanckCombined.fits"
    fn = f"{TENS_DIR}/{fn}"

    fh = fits.open(fn)
    ww = wcs.WCS(fh[0].header)

    return fh


def get_planck_data():
    center = coordinates.SkyCoord(0*u.deg, 0*u.deg, frame='galactic')

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
        # bandpass information here: https://wiki.cosmos.esa.int/planckpla2015/index.php/The_RIMO#HFI_2
        # adopt central frequency assuming spectral index = 3 giving nu_c for 100_avg = 104.225
        planck_image.header['REFFREQ'] = reffreq_hz
        planck_image.header['BMAJ'] = 9.65/60
        planck_image.header['BMIN'] = 9.65/60
        planck_image.header['BUNIT'] = 'K'
        planck_image.writeto(planckfn, overwrite=True)

    if not os.path.exists(planckfn_scaled):

        # scale planck data to match MGPS assuming alpha=3
        # mustang_reffreq = 92.44e9
        # revised: use alpha=0 instead
        planck_image.data = planck_image.data * (reffreq_hz / mustang_reffreq)**-3
        planck_image.header['REFFREQ'] = mustang_reffreq
        planck_image.writeto(planckfn_scaled, overwrite=True)


    fh[0].header['BMAJ'] = mustang_beam_fwhm.to(u.deg).value
    fh[0].header['BMIN'] = mustang_beam_fwhm.to(u.deg).value
    fh[0].header['REFFREQ'] = mustang_central_frequency.to(u.Hz).value
    fh[0].header['BUNIT'] = 'Jy/beam'

    # smooth 2024 MUSTANG data a tiny bit.  1 pixel = 2"
    fh[0].data = convolve_fft(fh[0].data, Gaussian2DKernel(1))

    # manually handle unit conversion - something in uvcombine broke?
    mustang_beam = Beam(mustang_beam_fwhm)
    jtok_mustang = mustang_beam.jtok(mustang_central_frequency)
    header_copy = fh[0].header.copy()
    header_copy['BUNIT'] = 'K'
    fh_K = fits.PrimaryHDU(fh[0].data*jtok_mustang, header=header_copy)

    rslt = feather_simple(fh_K, planckfn_scaled)

    fh[0].data = rslt.real
    fh[0].header['BUNIT'] = 'K'
    fh.writeto(mustang_plus_planck_fn, overwrite=True)