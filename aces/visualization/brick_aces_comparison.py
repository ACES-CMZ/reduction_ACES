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
brickfn = '/orange/adamginsburg/brick/alma/rathborne/brick.cont.alma.image.fits'
brickfeatherfn = '/orange/adamginsburg/brick/alma/rathborne/brick.cont.alma_sd.image.fits'
acesmosaicfn = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits'
acesmosaicfeatherfn = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'
mustangfn = '/orange/adamginsburg/ACES/TENS/SgrB2_flux_cut_dt6_filtp05to49_noShift_final_map_PlanckCombined.fits'

# Load data
brickfh = fits.open(brickfn)
brickfeatherfh = fits.open(brickfeatherfn)
acesmosaicfh = fits.open(acesmosaicfn)
acesmosaicfeatherfh = fits.open(acesmosaicfeatherfn)
brickcube = SpectralCube.read(brickfn)
mustangfh = fits.open(mustangfn)
mustang_reffreq = 87.85e9 * u.Hz # for alpha=0

# Reproject ACES to Brick
aces_to_brick, _ = reproject.reproject_interp(acesmosaicfh, brickcube[0].hdu.header) # Jy/beam
aces_to_brick_feather, _ = reproject.reproject_interp(acesmosaicfeatherfh, brickcube[0].hdu.header)

# Beam and convolution
BrickBeam = Beam.from_fits_header(brickfh[0].header)
ACESbeam = Beam.from_fits_header(acesmosaicfh[0].header)
ww = WCS(brickfh[0].header).celestial
pixscale = ww.proj_plane_pixel_area()**0.5
convkernel = ACESbeam.deconvolve(BrickBeam).as_kernel(pixscale)
brick_conv_ACES = convolve_fft(brickfh[0].data.squeeze(), convkernel)
brick_feather_conv_ACES = convolve_fft(brickfeatherfh[0].data.squeeze(), convkernel)

# Hide bad pixels
brick_conv_ACES[np.isnan(brickfh[0].data.squeeze())] = np.nan

# Conversion factors
ACES_reffreq = 97.21 * u.GHz
jtok_ACES = ACESbeam.jtok(ACES_reffreq)
nueff_brick = (87.2 + 89.1 + 99.1 + 101.1) / 4 * u.GHz
jtok_Brick = BrickBeam.jtok(nueff_brick)

# Generate figures
fig = pl.figure(figsize=(12, 11))
ax1 = pl.subplot(2, 2, 1, projection=ww)
brickdata = brickfh[0].data.squeeze() * jtok_Brick.value
norm = simple_norm(brickdata, stretch='linear', vmin=-0.01, vmax=0.055)
pl.imshow(brickdata[75:-75, 125:-125], cmap=mymap, norm=norm)

ax2 = pl.subplot(2, 2, 2, projection=ww)
pl.imshow(brick_conv_ACES[75:-75, 125:-125] * jtok_Brick.value, cmap=mymap, norm=norm)

ax3 = pl.subplot(2, 2, 3, projection=ww)
pl.imshow(aces_to_brick[75:-75, 125:-125] * jtok_ACES.value, cmap=mymap, norm=norm)

ax4 = pl.subplot(2, 2, 4, projection=ww)
diff = aces_to_brick * jtok_ACES.value - brick_conv_ACES * jtok_Brick.value
diff[np.isnan(brick_conv_ACES) | (brick_conv_ACES == 0) | np.isnan(brickdata)] = np.nan
pl.imshow(diff[75:-75, 125:-125], cmap='RdBu', norm=simple_norm(diff, vmin=-0.025, vmax=0.025, stretch='linear'))

# Add labels
ax1.text(0.95, 0.95, "Brick C0", transform=ax1.transAxes, horizontalalignment='right')
ax2.text(0.95, 0.95, "Brick C0 (smooth)", transform=ax2.transAxes, horizontalalignment='right')
ax3.text(0.95, 0.95, "ACES", transform=ax3.transAxes, horizontalalignment='right')
ax4.text(0.95, 0.95, "ACES - Bricksm", transform=ax4.transAxes, horizontalalignment='right')

# Format axes
cb1 = format_ax(ax1, label="", hidex=True)
cb3 = format_ax(ax3, label="")
cb2 = format_ax(ax2, label="T$_B$ [K]", cb=True, hidey=True, hidex=True)
cb4 = format_ax(ax4, label="ACES - C0", cb=True, hidey=True, hidex=False)

ax3.coords[0].set_ticklabel(rotation=15, pad=25)
ax4.coords[0].set_ticklabel(rotation=15, pad=25)

# Save figure
output_path = os.path.join(basepath, 'papers/continuum_data/figures/Brick_Cycle0_comparison_zoom.png')
pl.subplots_adjust(hspace=0.02, wspace=0.05)
pl.savefig(output_path, bbox_inches='tight', dpi=200)
print(f"Figure saved to {output_path}")

fig = pl.figure(figsize=(12, 11))
ax1 = pl.subplot(2, 2, 1, projection=ww[75:-75, 125:-125])
brickfeatherdata = brickfeatherfh[0].data.squeeze() * jtok_Brick.value
norm = simple_norm(brickfeatherdata, stretch='linear', vmin=-0.01, vmax=0.1) #max_percent=99.95, min_percent=1)
pl.imshow(brickfeatherdata[75:-75, 125:-125], cmap=mymap, norm=norm)
ax2 = pl.subplot(2, 2, 2, projection=ww[75:-75, 125:-125])
pl.imshow(brick_feather_conv_ACES[75:-75, 125:-125] * jtok_Brick.value, cmap=mymap, norm=norm)
ax3 = pl.subplot(2, 2, 3, projection=ww[75:-75, 125:-125])
pl.imshow(aces_to_brick_feather[75:-75, 125:-125], cmap=mymap, norm=norm)
ax4 = pl.subplot(2, 2, 4, projection=ww[75:-75, 125:-125])
diff = aces_to_brick_feather - brick_feather_conv_ACES*jtok_Brick.value
diff[np.isnan(brick_feather_conv_ACES) | (brick_feather_conv_ACES == 0) | np.isnan(brickfeatherdata)] = np.nan

brick_diff_limits = 0.06
pl.imshow(diff[75:-75, 125:-125], cmap='RdBu', norm=simple_norm(diff, vmin=-brick_diff_limits, vmax=brick_diff_limits, stretch='linear'))
#ax4.contour(diff, levels=[-5, -4, -3, -2, -1], colors=['w']*5, linewidths=[0.5]*5)
ax1.text(0.95, 0.95, "Feathered C0", transform=ax1.transAxes, horizontalalignment='right')
ax2.text(0.95, 0.95, "Feathered C0 (smooth)", transform=ax2.transAxes, horizontalalignment='right')
ax3.text(0.95, 0.95, "ACES", transform=ax3.transAxes, horizontalalignment='right')
ax4.text(0.95, 0.95, "ACES - smooth C0", transform=ax4.transAxes, horizontalalignment='right')

format_ax(ax1, label="", hidex=True)
format_ax(ax3, label="")
format_ax(ax2, label="T$_B$ [K]", cb=True, hidey=True, hidex=True)
format_ax(ax4, label="ACES - C0", cb=True, hidey=True, hidex=False)

# ignore = diff > -0.3
# diff[ignore] = np.nan
# ax4.imshow(diff, cmap='RdGy_r', norm=simple_norm(diff, vmin=-10, vmax=-0.3, stretch='linear'))

pl.subplots_adjust(hspace=0.02, wspace=0.05)
pl.savefig(f'{basepath}/papers/continuum_data/figures/BrickFeather_Cycle0_comparison_zoom.png', bbox_inches='tight', dpi=200)

pk = aces_to_brick_feather[75:-75, 125:-125].max()
ax4.contour(aces_to_brick_feather[75:-75, 125:-125], levels=np.array([0.1, 0.5, 0.9])*pk, colors=['g']*3, linewidths=[0.25]*3)
pl.savefig(f'{basepath}/papers/continuum_data/figures/BrickFeather_Cycle0_comparison_zoom_contours.png', bbox_inches='tight', dpi=200)

meerkatfh = fits.open('/orange/adamginsburg/cmz/meerkat/MEERKAT_StokesI_butreadablebyCARTA.fits')

target_hdu = brickcube[0].hdu.header
target_hdu['NAXIS1'] = 1600
target_hdu['CRPIX1'] = 801
target_hdu['NAXIS2'] = 1600
target_hdu['CRPIX2'] = 801
ww_big = WCS(target_hdu)

meerkat_to_brick, _ = reproject.reproject_interp(meerkatfh, target_hdu) # Jy/beam
aces_to_brick2, _ = reproject.reproject_interp(acesmosaicfeatherfh, target_hdu) # K
mustang_to_brick, _ = reproject.reproject_interp(mustangfh, target_hdu)

meerkat_beam = Beam.from_fits_header(meerkatfh[0].header)
meerkat_convkernel = meerkat_beam.deconvolve(ACESbeam).as_kernel(pixscale)
meerkat_reffreq = 1.4*u.GHz
meerkat_jtok = meerkat_beam.jtok(meerkat_reffreq)

mustang_beam = Beam.from_fits_header(mustangfh[0].header)
mustang_jtok = mustang_beam.jtok(mustang_reffreq) # MUSTANG is in K
meerkat_to_mustang_convkernel = mustang_beam.deconvolve(meerkat_beam).as_kernel(pixscale)

ACES_conv_MEER = convolve_fft(aces_to_brick2, meerkat_convkernel) # K
meerkat_conv_MUSTANG = convolve_fft(meerkat_to_brick, meerkat_to_mustang_convkernel)

effective_index = -0.1
frequency_scale_meeralma = (ACES_reffreq / meerkat_reffreq)**(2+effective_index)

mustang_meerkat_ratio = mustang_to_brick / (meerkat_conv_MUSTANG * meerkat_jtok.value * frequency_scale_meeralma)
mustang_meerkat_difference = mustang_to_brick - (meerkat_conv_MUSTANG * meerkat_jtok.value * frequency_scale_meeralma)

alpha_mustang_meerkat = np.log(mustang_meerkat_ratio) / np.log(mustang_reffreq / meerkat_reffreq)

# ACES_conv_MEER is in Kelvin but is now in a MEERKAT beam, so we need to go K->Jy in Meerkat beam units
aces_meerkatbeam_jtok = meerkat_beam.jtok(ACES_reffreq)
alpha_aces_meerkat = np.log(ACES_conv_MEER / aces_meerkatbeam_jtok.value / (meerkat_to_brick)) / np.log(ACES_reffreq / meerkat_reffreq)


fig = pl.figure(figsize=(12, 11))
slc = slice(475, -475), slice(500, -500)
slc = slice(0, None), slice(0, None)
ax1 = pl.subplot(2, 2, 1, projection=ww_big[slc])
norm = simple_norm(meerkat_to_brick * meerkat_jtok.value / frequency_scale_meeralma,
                   stretch='linear', vmin=-0.01, vmax=0.1) #max_percent=99.95, min_percent=1)
pl.imshow((meerkat_to_brick * meerkat_jtok.value / frequency_scale_meeralma)[slc],
          cmap=mymap, norm=norm)
ax2 = pl.subplot(2, 2, 2, projection=ww_big[slc])
pl.imshow(ACES_conv_MEER[slc], cmap=mymap, norm=norm)

ax3 = pl.subplot(2, 2, 4, projection=ww_big[slc])
pl.imshow(alpha_aces_meerkat[slc], cmap='Spectral', norm=simple_norm(alpha_aces_meerkat[slc], vmin=-0.3, vmax=0.8))

ax4 = pl.subplot(2, 2, 3, projection=ww_big[slc])
diff = -(meerkat_to_brick * meerkat_jtok.value / frequency_scale_meeralma
         - ACES_conv_MEER)
pl.imshow(diff[slc], cmap=mymap, norm=norm)

ax1.text(0.95, 0.95, f"MEERKAT", transform=ax1.transAxes, horizontalalignment='right')
ax2.text(0.95, 0.95, "ACES smooth", transform=ax2.transAxes, horizontalalignment='right',)
ax3.text(0.95, 0.95, "$\\alpha$(ACESsm / MEERKAT)", transform=ax3.transAxes, horizontalalignment='right')
ax4.text(0.95, 0.95, f"ACESsm-MEERKAT", transform=ax4.transAxes, horizontalalignment='right')

format_ax(ax1, label="", cb=False, hidex=True)
format_ax(ax4, label="", cb=False)
format_ax(ax2, label="T$_B$ [K]", cb=True, hidey=True, hidex=True)
format_ax(ax3, label="$\\alpha$", cb=True, hidey=True, hidex=False)

ax3.coords[0].set_ticklabel(rotation=25, pad=35)
ax4.coords[0].set_ticklabel(rotation=25, pad=35)

pl.subplots_adjust(hspace=0.02, wspace=0.02)
pl.savefig(f'{basepath}/papers/continuum_data/figures/BrickMEERKAT_comparison_large.png',
           bbox_inches='tight', dpi=200)
print(f'Effective spectral index assumed={effective_index:0.3f}  frequency scaling factor={frequency_scale_meeralma:0.3f}')


convkernel_bricktomeer = meerkat_beam.deconvolve(BrickBeam).as_kernel(pixscale)
brick_feather_conv_MEER = convolve_fft(brickfeatherfh[0].data.squeeze(),
                                       convkernel_bricktomeer)
brick_feather_conv_MEER[np.isnan(brickfh[0].data.squeeze())] = np.nan

fig = pl.figure(figsize=(10, 17))
slc = slice(475, -475), slice(525, -525)
ax1 = pl.subplot(3, 2, 1, projection=ww_big[slc])
norm = simple_norm(meerkat_to_brick * meerkat_jtok.value / frequency_scale_meeralma,
                   stretch='linear', vmin=-0.01, vmax=0.1) #max_percent=99.95, min_percent=1)
pl.imshow((meerkat_to_brick * meerkat_jtok.value / frequency_scale_meeralma)[slc],
          cmap=mymap, norm=norm)
ax2 = pl.subplot(3, 2, 2, projection=ww_big[slc])
pl.imshow(ACES_conv_MEER[slc], cmap=mymap, norm=norm)
ax3 = pl.subplot(3, 2, 3, projection=ww_big[slc])
diff = -(meerkat_to_brick * meerkat_jtok.value / frequency_scale_meeralma
         - ACES_conv_MEER)
pl.imshow(diff[slc], cmap=mymap, norm=norm)

ax4 = pl.subplot(3, 2, 4, projection=ww_big[slc])
#pl.imshow(aces_to_brick2[slc], cmap=mymap, norm=norm)
pl.imshow(alpha_aces_meerkat[slc], cmap='Spectral', norm=simple_norm(alpha_aces_meerkat[slc], vmin=-0.2, vmax=0.7))

ax6 = pl.subplot(3, 2, 5, projection=ww[75:-75, 125:-125])
diff = -(meerkat_to_brick * meerkat_jtok.value / frequency_scale_meeralma
         - ACES_conv_MEER
        )[400:-400, 400:-400] - brick_feather_conv_MEER*jtok_Brick.value
pl.imshow(diff[75:-75, 125:-125], cmap='RdBu',
          norm=simple_norm(diff, vmin=-brick_diff_limits, vmax=brick_diff_limits, stretch='linear'))

ax5 = pl.subplot(3, 2, 6, projection=ww[75:-75, 125:-125])
pl.imshow(brick_feather_conv_MEER[75:-75, 125:-125] * jtok_Brick.value,
          cmap=mymap, norm=norm)

ax1.text(0.95, 0.95, f"MEERKAT", transform=ax1.transAxes, horizontalalignment='right')
ax2.text(0.95, 0.95, "ACES smooth", transform=ax2.transAxes, horizontalalignment='right', color='w')
ax3.text(0.95, 0.95, f"ACESsm-MEERKAT", transform=ax3.transAxes, horizontalalignment='right', color='w')
ax4.text(0.95, 0.95, "$\\alpha$ ACES/MEERKAT", transform=ax4.transAxes, horizontalalignment='right', color='k')
ax5.text(0.95, 0.95, "BrickC0", transform=ax5.transAxes, horizontalalignment='right')
ax6.text(0.95, 0.95, f"ACESsm-MEERKAT-C0", transform=ax6.transAxes, horizontalalignment='right')

format_ax(ax1, label="", cb=False, hidex=True)
format_ax(ax3, label="", cb=False, hidex=True)
format_ax(ax2, label="T$_B$ [K]", cb=True, hidey=True, hidex=True)
format_ax(ax4, label="$\\alpha$", cb=True, hidey=True, hidex=True)
format_ax(ax5, label="T$_B$ [K]", cb=True, hidey=True, hidex=False)
format_ax(ax6, label="", cb=False, hidey=False, hidex=False)

ax5.coords[0].set_ticklabel(rotation=25, pad=35)
ax6.coords[0].set_ticklabel(rotation=25, pad=35)

pl.tight_layout()
pl.subplots_adjust(hspace=0.02, wspace=0.0)
pl.savefig(f'{basepath}/papers/continuum_data/figures/BrickMEERKAT_comparison_large_withC0.png',
           bbox_inches='tight', dpi=200)
print(f'Effective spectral index assumed={effective_index:0.3f}')





fig = pl.figure(figsize=(12, 11))
ax1 = pl.subplot(2, 2, 1, projection=ww_big)
norm = simple_norm(meerkat_conv_MUSTANG * meerkat_jtok.value * frequency_scale_meeralma,
                   stretch='linear', vmin=-10, vmax=50) #max_percent=99.95, min_percent=1)
pl.imshow(meerkat_conv_MUSTANG * meerkat_jtok.value * frequency_scale_meeralma,
          cmap=mymap, norm=norm)
ax2 = pl.subplot(2, 2, 2, projection=ww_big)
pl.imshow(mustang_to_brick, cmap=mymap, norm=norm) #simple_norm(mustang_to_brick, stretch='linear', max_percent=99.5, min_percent=1))
ax3 = pl.subplot(2, 2, 3, projection=ww_big)
pl.imshow(mustang_meerkat_difference, 
          norm=simple_norm(mustang_meerkat_difference, vmin=-5, vmax=5),
          cmap=mymap)
          #norm=simple_norm(mustang_meerkat_ratio, stretch='linear', vmin=-0.001, vmax=0.005))
ax4 = pl.subplot(2, 2, 4, projection=ww_big)
pl.imshow(alpha_mustang_meerkat, norm=simple_norm(alpha_mustang_meerkat, stretch='linear', vmin=-1.25, vmax=0.1))

ax1.text(0.95, 0.95, "MEERKAT", transform=ax1.transAxes, horizontalalignment='right')
ax2.text(0.95, 0.95, "MUSTANG", transform=ax2.transAxes, horizontalalignment='right')
ax3.text(0.95, 0.95, "MUSTANG - MEERKAT", transform=ax3.transAxes, horizontalalignment='right')
ax4.text(0.95, 0.95, "$\\alpha_{MUSTANG/MEERKAT}$", transform=ax4.transAxes, horizontalalignment='right')

format_ax(ax1, label="", cb=False, hidex=True)
format_ax(ax3, label="", cb=False)
format_ax(ax2, label="T$_B$ [K]", cb=True, hidey=True, hidex=True)
format_ax(ax4, label="$\\alpha$", cb=True, hidey=True, hidex=False)

ax3.coords[0].set_ticklabel(rotation=25, pad=35)
ax4.coords[0].set_ticklabel(rotation=25, pad=35)

pl.subplots_adjust(hspace=0.02, wspace=0.02)
pl.savefig(f'{basepath}/papers/continuum_data/figures/BrickMEERKATvsMUSTANG_comparison.png',
           bbox_inches='tight', dpi=200)