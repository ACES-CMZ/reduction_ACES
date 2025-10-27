from spectral_cube import SpectralCube
from astropy.io import fits
import pylab as pl
from radio_beam import Beam
from astropy.convolution import convolve, convolve_fft
from astropy.wcs import WCS
from aces import conf

import reproject

import matplotlib.colors as mcolors
import numpy as np
from astropy.visualization import simple_norm
from astropy import units as u

from mpl_toolkits.axes_grid1 import make_axes_locatable

basepath = conf.basepath

pl.rcParams['font.size'] = 16

def format_ax(ax, label="", cb=True, hidey=False, hidex=False, cbar_size="5%"):
    ax.coords[0].set_axislabel("Right Ascension")
    ax.coords[1].set_axislabel("Declination")
    im = ax.images[0]

    if cb:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size=cbar_size, pad=0.05)
        cax.set_xticks([])
        cax.set_xlabel("")
        cax.set_xticklabels([])

        cb = pl.colorbar(im, cax=cax)
        cax.coords[0].set_axislabel('')
        cax.coords[0].set_ticklabel_visible(False)
        cax.coords[0].set_ticks_visible(False)
        cax.coords[1].set_ticks_position('r')
        cax.coords[1].set_ticklabel_position('r')
        if label:
            cax.coords[1].set_axislabel_position('r')
        else:
            cax.coords[1].set_axislabel_position('')
        cax.coords[1].set_axislabel(" ")
        cax.coords[1].set_axislabel(label)
    if hidey:
        ax.coords[1].set_ticklabel_visible(False)
        ax.coords[1].set_axislabel("")
    if hidex:
        ax.coords[0].set_ticklabel_visible(False)
        ax.coords[0].set_axislabel("")

    return cb


def compare_to_other_data(otherfn,
                          name,
                          nueff_other,
                          acesfn=f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits',
                          norm_kwargs={},
                          hspace=0.01,
                          wspace=0.1,
                          figsize=(12, 10),
                          diffnorm_kwargs={'min_percent': 0.1, 'max_percent': 99.9, 'stretch': 'linear'},
                          ticklabel_pad=None,
                          ticklabel_rotation=None,
                          ):

    # colormap setup
    colors1 = pl.cm.gray_r(np.linspace(0., 1, 128))
    colors2 = pl.cm.hot(np.linspace(0, 1, 128))

    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    otherfh = fits.open(otherfn)
    #othercube = SpectralCube.read(otherfn)

    acesmosaicfh = fits.open(acesfn)

    aces_to_other, _ = reproject.reproject_interp(acesmosaicfh, otherfh[0].header)

    otherBeam = Beam.from_fits_header(otherfh[0].header)
    ACESbeam = Beam.from_fits_header(acesmosaicfh[0].header)
    ww = WCS(otherfh[0].header)
    pixscale = ww.proj_plane_pixel_area()**0.5
    convkernel = ACESbeam.deconvolve(otherBeam).as_kernel(pixscale)
    other_conv_ACES = convolve_fft(otherfh[0].data.squeeze(), convkernel)

    jtok_ACES = ACESbeam.jtok(97.21 * u.GHz)
    jtok_other = otherBeam.jtok(nueff_other)

    beam_ratio = (ACESbeam.sr / otherBeam.sr)  # noqa [I want beam_ratio defined]
    fig = pl.figure(figsize=figsize)  # noqa [I want fig defined]
    ax1 = pl.subplot(2, 2, 1, projection=ww)
    otherdata = otherfh[0].data.squeeze() * jtok_other.value
    norm = simple_norm(otherdata, **norm_kwargs)
    pl.imshow(otherdata, cmap=mymap, norm=norm)
    ax2 = pl.subplot(2, 2, 2, projection=ww)
    pl.imshow(other_conv_ACES * jtok_other.value, cmap=mymap, norm=norm)
    ax3 = pl.subplot(2, 2, 3, projection=ww)
    pl.imshow(aces_to_other * jtok_ACES.value, cmap=mymap, norm=norm)
    ax4 = pl.subplot(2, 2, 4, projection=ww)
    diff = aces_to_other * jtok_ACES.value - other_conv_ACES * jtok_other.value
    diff[np.isnan(other_conv_ACES) | (other_conv_ACES == 0) | np.isnan(otherdata)] = np.nan
    pl.imshow(diff, cmap='RdBu', norm=simple_norm(diff, **diffnorm_kwargs))

    format_ax(ax1, label="", hidex=True)
    format_ax(ax3, label="")
    format_ax(ax2, label="T$_B$ [K]", cb=True, hidey=True, hidex=True)
    format_ax(ax4, label=f"ACES - {name}", cb=True, hidey=True, hidex=False)

    ax1.text(0.95, 0.95, f"{name}", transform=ax1.transAxes, horizontalalignment='right')
    ax2.text(0.95, 0.95, f"{name} (smooth)", transform=ax2.transAxes, horizontalalignment='right')
    ax3.text(0.95, 0.95, "ACES", transform=ax3.transAxes, horizontalalignment='right')
    ax4.text(0.95, 0.95, f"ACES - {name}sm", transform=ax4.transAxes, horizontalalignment='right')

    for lon in (ax3.coords[0], ax4.coords[0]):
        if ticklabel_rotation is not None:
            lon.set_ticklabel(rotation=ticklabel_rotation, pad=ticklabel_pad)

    pl.subplots_adjust(hspace=hspace, wspace=wspace)
    pl.savefig(f'{basepath}/papers/continuum_data/figures/{name}_comparison.png', bbox_inches='tight', dpi=200)