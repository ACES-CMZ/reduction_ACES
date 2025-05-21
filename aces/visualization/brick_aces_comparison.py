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

# Paths and filenames
basepath = conf.basepath
brickfn = '/orange/adamginsburg/brick/alma/rathborne/brick.cont.alma.image.fits'
brickfeatherfn = '/orange/adamginsburg/brick/alma/rathborne/brick.cont.alma_sd.image.fits'
acesmosaicfn = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits'
acesmosaicfeatherfn = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'

# Load data
brickfh = fits.open(brickfn)
brickfeatherfh = fits.open(brickfeatherfn)
acesmosaicfh = fits.open(acesmosaicfn)
acesmosaicfeatherfh = fits.open(acesmosaicfeatherfn)
brickcube = SpectralCube.read(brickfn)

# Reproject ACES to Brick
aces_to_brick, _ = reproject.reproject_interp(acesmosaicfh, brickcube[0].hdu.header)
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
jtok_ACES = ACESbeam.jtok(97.21 * u.GHz)
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
pl.imshow(diff[75:-75, 125:-125], cmap='RdBu', norm=simple_norm(diff, vmin=-0.04, vmax=0.04, stretch='linear'))
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

meerkat_to_brick, _ = reproject.reproject_interp(meerkatfh, target_hdu)
aces_to_brick2, _ = reproject.reproject_interp(acesmosaicfeatherfh, target_hdu)

meerkat_beam = Beam.from_fits_header(meerkatfh[0].header)
meerkat_convkernel = meerkat_beam.deconvolve(ACESbeam).as_kernel(pixscale)
meerkat_jtok = meerkat_beam.jtok(1.4*u.GHz)

ACES_conv_MEER = convolve_fft(aces_to_brick2, meerkat_convkernel)

fig = pl.figure(figsize=(12, 11))
ax1 = pl.subplot(2, 2, 1, projection=ww_big)
arbitrary_scale = 0.00015 #(95/1.4)**-0.8
effective_index = np.log(arbitrary_scale) / np.log(95/1.4)
norm = simple_norm(meerkat_to_brick * meerkat_jtok.value * arbitrary_scale,
                   stretch='linear', vmin=-0.01, vmax=0.1) #max_percent=99.95, min_percent=1)
pl.imshow(meerkat_to_brick * meerkat_jtok.value * arbitrary_scale,
          cmap=mymap, norm=norm)
ax2 = pl.subplot(2, 2, 2, projection=ww_big)
pl.imshow(ACES_conv_MEER, cmap=mymap, norm=norm)
ax3 = pl.subplot(2, 2, 3, projection=ww_big)
pl.imshow(aces_to_brick2, cmap=mymap, norm=norm)
ax4 = pl.subplot(2, 2, 4, projection=ww_big)
diff = -(meerkat_to_brick * meerkat_jtok.value * arbitrary_scale
         - ACES_conv_MEER)
pl.imshow(diff, cmap=mymap, norm=norm)

ax1.text(0.95, 0.95, f"MEERKAT*{arbitrary_scale}", transform=ax1.transAxes, horizontalalignment='right')
ax2.text(0.95, 0.95, "ACES smooth", transform=ax2.transAxes, horizontalalignment='right')
ax3.text(0.95, 0.95, "ACES", transform=ax3.transAxes, horizontalalignment='right')
ax4.text(0.95, 0.95, f"ACES-MEERKAT*{arbitrary_scale}", transform=ax4.transAxes, horizontalalignment='right')

format_ax(ax1, label="", cb=False, hidex=True)
format_ax(ax3, label="", cb=False)
format_ax(ax2, label="T$_B$ [K]", cb=True, hidey=True, hidex=True)
format_ax(ax4, label="T$_B$ [K]", cb=True, hidey=True, hidex=False)

ax3.coords[0].set_ticklabel(rotation=25, pad=30)
ax4.coords[0].set_ticklabel(rotation=25, pad=30)

pl.subplots_adjust(hspace=0.02, wspace=0.02)
pl.savefig(f'{basepath}/papers/continuum_data/figures/BrickMEERKAT_comparison_large.png',
           bbox_inches='tight', dpi=200)
print(f'Effective spectral index assumed={effective_index:0.3f}')


convkernel_bricktomeer = meerkat_beam.deconvolve(BrickBeam).as_kernel(pixscale)
brick_feather_conv_MEER = convolve_fft(brickfeatherfh[0].data.squeeze(),
                                       convkernel_bricktomeer)
brick_feather_conv_MEER[np.isnan(brickfh[0].data.squeeze())] = np.nan

fig = pl.figure(figsize=(12, 17))
slc = slice(475, -475), slice(525, -525)
ax1 = pl.subplot(3, 2, 1, projection=ww_big[slc])
arbitrary_scale = 0.00015 #(95/1.4)**-0.8
effective_index = np.log(arbitrary_scale) / np.log(95/1.4)
norm = simple_norm(meerkat_to_brick * meerkat_jtok.value * arbitrary_scale,
                   stretch='linear', vmin=-0.01, vmax=0.1) #max_percent=99.95, min_percent=1)
pl.imshow((meerkat_to_brick * meerkat_jtok.value * arbitrary_scale)[slc],
          cmap=mymap, norm=norm)
ax2 = pl.subplot(3, 2, 2, projection=ww_big[slc])
pl.imshow(ACES_conv_MEER[slc], cmap=mymap, norm=norm)
ax3 = pl.subplot(3, 2, 3, projection=ww_big[slc])
pl.imshow(aces_to_brick2[slc], cmap=mymap, norm=norm)
ax4 = pl.subplot(3, 2, 4, projection=ww_big[slc])
diff = -(meerkat_to_brick * meerkat_jtok.value * arbitrary_scale
         - ACES_conv_MEER)
pl.imshow(diff[slc], cmap=mymap, norm=norm)

ax5 = pl.subplot(3, 2, 5, projection=ww[75:-75, 125:-125])
diff = -(meerkat_to_brick * meerkat_jtok.value * arbitrary_scale
         - ACES_conv_MEER
        )[400:-400, 400:-400] - brick_feather_conv_MEER*jtok_Brick.value
pl.imshow(diff[75:-75, 125:-125], cmap='RdBu',
          norm=simple_norm(diff, vmin=-0.04, vmax=0.04, stretch='linear'))

ax6 = pl.subplot(3, 2, 6, projection=ww[75:-75, 125:-125])
pl.imshow(brick_feather_conv_MEER[75:-75, 125:-125] * jtok_Brick.value,
          cmap=mymap, norm=norm)

ax1.text(0.95, 0.95, f"MEERKAT*{arbitrary_scale}", transform=ax1.transAxes, horizontalalignment='right')
ax2.text(0.95, 0.95, "ACES smooth", transform=ax2.transAxes, horizontalalignment='right')
ax3.text(0.95, 0.95, "ACES", transform=ax3.transAxes, horizontalalignment='right')
ax4.text(0.95, 0.95, f"ACES-MEERKAT*{arbitrary_scale}", transform=ax4.transAxes, horizontalalignment='right')
ax5.text(0.95, 0.95, f"ACES-MEERKAT*{arbitrary_scale}-Brick", transform=ax5.transAxes, horizontalalignment='right')
ax6.text(0.95, 0.95, "BrickC0", transform=ax6.transAxes, horizontalalignment='right')

format_ax(ax1, label="", cb=False, hidex=True)
format_ax(ax3, label="", cb=False, hidex=True)
format_ax(ax2, label="T$_B$ [K]", cb=True, hidey=True, hidex=True)
format_ax(ax4, label="T$_B$ [K]", cb=True, hidey=True, hidex=True)
format_ax(ax5, label="", cb=True, hidey=False, hidex=False)
format_ax(ax6, label="T$_B$ [K]", cb=True, hidey=True, hidex=False)

#ax5.coords[0].set_ticklabel(rotation=25, pad=30)
#ax4.coords[0].set_ticklabel(rotation=25, pad=30)

pl.subplots_adjust(hspace=0.02, wspace=0.04)
pl.savefig(f'{basepath}/papers/continuum_data/figures/BrickMEERKAT_comparison_large_withC0.png',
           bbox_inches='tight', dpi=200)
print(f'Effective spectral index assumed={effective_index:0.3f}')
