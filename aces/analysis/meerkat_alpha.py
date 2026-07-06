import os
from spectral_cube import SpectralCube
from astropy.io import fits
import pylab as pl
from radio_beam import Beam
from astropy.convolution import convolve_fft
from astropy.wcs import WCS
from aces import conf
from aces.visualization.figure_configuration import (
    mymap, format_ax
)
import reproject
import numpy as np
from astropy.visualization import simple_norm
from astropy import units as u
from aces.visualization.merge_mustang import get_mustang_data

# Paths and filenames
basepath = conf.basepath
acesmosaicfn = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits'
acesmosaicfeatherfn = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'
mustangfn = '/orange/adamginsburg/ACES/TENS/SgrB2_flux_cut_dt6_filtp05to49_noShift_final_map_PlanckCombined.fits'
meerkatfh = fits.open('/orange/adamginsburg/cmz/meerkat/MEERKAT_StokesI_butreadablebyCARTA.fits')

# Load data
acesmosaicfh = fits.open(acesmosaicfn)
acesmosaicfeatherfh = fits.open(acesmosaicfeatherfn)
mustangfh = fits.open(mustangfn)
mustang_reffreq = 87.85e9 * u.Hz # for alpha=0


# Beam and convolution
ACESbeam = Beam.from_fits_header(acesmosaicfh[0].header)
ww = WCS(acesmosaicfeatherfh[0].header).celestial
pixscale = ww.proj_plane_pixel_area()**0.5

# Conversion factors
ACES_reffreq = 97.21 * u.GHz
jtok_ACES = ACESbeam.jtok(ACES_reffreq)

target_hdu = acesmosaicfh[0].header
# target_hdu['NAXIS1'] = 1600
# target_hdu['CRPIX1'] = 801
# target_hdu['NAXIS2'] = 1600
# target_hdu['CRPIX2'] = 801
ww_big = WCS(target_hdu)


# Reproject MEERKAT to ACES.  MeerKAT FITS data are Jy/beam; attach units now.
meerkat_to_aces_arr, _ = reproject.reproject_interp(meerkatfh, target_hdu)
meerkat_to_aces = u.Quantity(meerkat_to_aces_arr, u.Jy / u.beam)
meerkat_beam = Beam.from_fits_header(meerkatfh[0].header)
meerkat_convkernel = meerkat_beam.deconvolve(ACESbeam).as_kernel(pixscale)
meerkat_reffreq = 1.284*u.GHz  # SARAO MeerKAT galactic centre mosaic (Heywood+ 2022) L-band center
meerkat_jtok = meerkat_beam.jtok(meerkat_reffreq)

# not used mustang_beam = Beam.from_fits_header(mustangfh[0].header)
# not used mustang_jtok = mustang_beam.jtok(mustang_reffreq) # MUSTANG is in K
# not used meerkat_to_mustang_convkernel = mustang_beam.deconvolve(meerkat_beam).as_kernel(pixscale)

# ACES feathered mosaic has BUNIT='K' (see joint_deconvolution/feather_continuum.py).
# convolve_fft strips units; we re-attach K immediately so downstream math stays unit-aware.
ACES_conv_MEER = u.Quantity(
    convolve_fft(acesmosaicfeatherfh[0].data, meerkat_convkernel),
    u.K,
)
# not used meerkat_conv_MUSTANG = convolve_fft(meerkat_to_aces, meerkat_to_mustang_convkernel)

effective_index = -0.1
frequency_scale_meeralma = (ACES_reffreq / meerkat_reffreq)**(2+effective_index)

# not used mustang_meerkat_ratio = mustang_to_aces / (meerkat_conv_MUSTANG * meerkat_jtok.value * frequency_scale_meeralma)
# not used mustang_meerkat_difference = mustang_to_aces - (meerkat_conv_MUSTANG * meerkat_jtok.value * frequency_scale_meeralma)

# alpha_mustang_meerkat = np.log(mustang_meerkat_ratio) / np.log(mustang_reffreq / meerkat_reffreq)

# ACES_conv_MEER is in K at the ACES freq, evaluated in a MEERKAT-sized beam.
# Divide by the MEERKAT-beam Jy<->K factor at the ACES frequency to get Jy/beam
# (in the MEERKAT beam) at the ACES freq, then take the dimensionless flux ratio
# vs the MEERKAT image (Jy/beam at MEERKAT freq) to compute the spectral index.
aces_meerkatbeam_jtok = meerkat_beam.jtok(ACES_reffreq)  # K per (Jy/beam)
aces_jybeam_at_meerkatbeam = (ACES_conv_MEER / aces_meerkatbeam_jtok).to(u.Jy / u.beam)
flux_ratio = (aces_jybeam_at_meerkatbeam / meerkat_to_aces).decompose()
# Guard: ratio must be truly dimensionless before log
assert flux_ratio.unit == u.dimensionless_unscaled, f"non-dimensionless ratio: {flux_ratio.unit}"
log_freq_ratio = np.log((ACES_reffreq / meerkat_reffreq).decompose().value)
alpha_aces_meerkat = np.log(flux_ratio.value) / log_freq_ratio

hdu = fits.PrimaryHDU(data=alpha_aces_meerkat, header=target_hdu)
hdu.writeto(f'{basepath}/mosaics/continuum/alpha_aces_meerkat.fits', overwrite=True)