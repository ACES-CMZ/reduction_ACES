import os
from spectral_cube import SpectralCube
from astropy.io import fits
import pylab as pl
from radio_beam import Beam
from astropy.convolution import convolve_fft
from astropy.wcs import WCS
from aces import conf
from aces.analysis.figure_scripts.figure_configuration import (
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
