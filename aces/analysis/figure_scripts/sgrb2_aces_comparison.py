#!/usr/bin/env python
"""
Script to compare ACES and SgrB2 2013 continuum data, creating comparison plots
between the two datasets.
"""

import numpy as np
import pylab as pl
import matplotlib.colors as mcolors
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.visualization import simple_norm
from astropy.convolution import convolve, convolve_fft
from radio_beam import Beam
import reproject

from aces.analysis.figure_scripts import figure_configuration
from aces.analysis.figure_scripts.figure_configuration import mymap, mymap2, distance

from aces import conf

basepath = conf.basepath


def load_data():
    """Load the SgrB2 and ACES mosaic data"""
    sgrb2fh = fits.open('/orange/adamginsburg/sgrb2/2013.1.00269.S/continuum/SgrB2_selfcal_full_TCTE7m_try2_selfcal6_ampphase_deeper_mask1.5mJy.image.tt0.pbcor.fits')
    acesmosaicfh = fits.open('/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits')
    acesmosaicfeatherfh = fits.open('/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits')
    
    return sgrb2fh, acesmosaicfh, acesmosaicfeatherfh


def reproject_data(sgrb2fh, acesmosaicfh, acesmosaicfeatherfh):
    """Reproject ACES data to the SgrB2 grid"""
    aces_to_sgr, _ = reproject.reproject_interp(acesmosaicfh, sgrb2fh[0].header)
    aces_to_sgr_feather, _ = reproject.reproject_interp(acesmosaicfeatherfh, sgrb2fh[0].header)
    
    return aces_to_sgr, aces_to_sgr_feather


def convolve_data(sgrb2fh, acesmosaicfh):
    """Convolve SgrB2 data to ACES resolution"""
    SgrB2Beam = Beam.from_fits_header(sgrb2fh[0].header)
    ACESbeam = Beam.from_fits_header(acesmosaicfh[0].header)
    ww = WCS(sgrb2fh[0].header)
    pixscale = ww.proj_plane_pixel_area()**0.5
    convkernel = ACESbeam.deconvolve(SgrB2Beam).as_kernel(pixscale)
    sgr_conv_ACES = convolve_fft(sgrb2fh[0].data, convkernel)
    
    return sgr_conv_ACES, SgrB2Beam, ACESbeam, ww


def calculate_jtok(ACESbeam, SgrB2Beam):
    """Calculate Jy/beam to K conversion factors"""
    jtok_ACES = ACESbeam.jtok(97.21*u.GHz)
    nueff_sgrb2 = (89.48 + 91.28 + 101.37 + 103.23)/4*u.GHz
    jtok_SgrB2 = SgrB2Beam.jtok(nueff_sgrb2)
    
    return jtok_ACES, jtok_SgrB2


def plot_comparison_zoom(sgrb2fh, sgr_conv_ACES, aces_to_sgr, ww, jtok_SgrB2, jtok_ACES, mymap, feathered=False):
    """Create zoomed comparison figure"""
    fig = pl.figure(figsize=(12, 11))
    
    # First panel: SgrB2 2013 data
    ax1 = pl.subplot(2, 2, 1, projection=ww)
    sgrb2data = sgrb2fh[0].data * jtok_SgrB2.value
    norm = simple_norm(sgrb2data, stretch='asinh', max_percent=99.95, min_percent=1)
    ax1.imshow(sgrb2data[1024:-1024, 1024:-1024], cmap=mymap, norm=norm)
    ax1.contour(sgrb2data[1024:-1024, 1024:-1024], colors=['k']*5, levels=[10, 15, 20])
    
    # Second panel: Smoothed SgrB2 data
    ax2 = pl.subplot(2, 2, 2, projection=ww)
    ax2.imshow(sgr_conv_ACES[1024:-1024, 1024:-1024] * jtok_SgrB2.value, cmap=mymap, norm=norm)
    ax2.contour(sgr_conv_ACES[1024:-1024, 1024:-1024] * jtok_SgrB2.value, colors=['k']*5, levels=[10, 15, 20])
    
    # Third panel: ACES data
    ax3 = pl.subplot(2, 2, 3, projection=ww)
    ax3.imshow(aces_to_sgr[1024:-1024, 1024:-1024] * jtok_ACES.value, cmap=mymap, norm=norm)
    ax3.contour(aces_to_sgr[1024:-1024, 1024:-1024] * jtok_ACES.value, colors=['k']*5, levels=[10, 15, 20])
    
    # Fourth panel: Difference
    ax4 = pl.subplot(2, 2, 4, projection=ww)
    diff = aces_to_sgr*jtok_ACES.value - sgr_conv_ACES*jtok_SgrB2.value
    pl.imshow(diff[1024:-1024, 1024:-1024], cmap='RdBu', norm=simple_norm(diff, vmin=-0.5, vmax=0.5, stretch='linear'))
    ax4.contour(diff[1024:-1024, 1024:-1024], levels=[-5, -4, -3, -2, -1], colors=['w']*5, linewidths=[0.5]*5)
    
    # Add labels
    ax1.text(0.95, 0.95, "Sgr B2 2013", transform=ax1.transAxes, horizontalalignment='right')
    ax2.text(0.95, 0.95, "Sgr B2 2013 (smooth)", transform=ax2.transAxes, horizontalalignment='right')
    ax3.text(0.95, 0.95, "ACES", transform=ax3.transAxes, horizontalalignment='right')
    ax4.text(0.95, 0.95, "ACES - SgrB2sm", transform=ax4.transAxes, horizontalalignment='right')
    
    # Format axes
    figure_configuration.format_ax(ax1, label="", hidex=True)
    figure_configuration.format_ax(ax3, label="")
    figure_configuration.format_ax(ax2, label="T$_B$ [K]", cb=True, hidey=True, hidex=True)
    figure_configuration.format_ax(ax4, label="ACES - 2013", cb=True, hidey=True, hidex=False)
    
    pl.subplots_adjust(hspace=0.01, wspace=0.05)
    
    # Save figure
    suffix = "feather_" if feathered else ""
    pl.savefig(f'{basepath}/papers/continuum_data/figures/SgrB2_2013_{suffix}comparison_zoom.png', 
               bbox_inches='tight', dpi=200)


def plot_comparison_full(sgrb2fh, sgr_conv_ACES, aces_to_sgr, ww, jtok_SgrB2, jtok_ACES, mymap):
    """Create full (non-zoomed) comparison figure"""
    fig = pl.figure(figsize=(12, 10))
    
    # First panel: SgrB2 2013 data
    ax1 = pl.subplot(2, 2, 1, projection=ww)
    sgrb2data = sgrb2fh[0].data * jtok_SgrB2.value
    norm = simple_norm(sgrb2data, stretch='asinh', max_percent=99.95, min_percent=1)
    pl.imshow(sgrb2data, cmap=mymap, norm=norm)
    
    # Second panel: Smoothed SgrB2 data
    ax2 = pl.subplot(2, 2, 2, projection=ww)
    pl.imshow(sgr_conv_ACES * jtok_SgrB2.value, cmap=mymap, norm=norm)
    
    # Third panel: ACES data
    ax3 = pl.subplot(2, 2, 3, projection=ww)
    pl.imshow(aces_to_sgr * jtok_ACES.value, cmap=mymap, norm=norm)
    
    # Fourth panel: Difference
    ax4 = pl.subplot(2, 2, 4, projection=ww)
    diff = aces_to_sgr*jtok_ACES.value - sgr_conv_ACES*jtok_SgrB2.value
    pl.imshow(diff, cmap='RdBu', norm=simple_norm(diff, min_percent=0.1, max_percent=99.0, stretch='asinh'))
    
    # Format axes
    figure_configuration.format_ax(ax1, label="")
    figure_configuration.format_ax(ax3, label="")
    figure_configuration.format_ax(ax2, label="T$_B$ [K]", cb=True)
    figure_configuration.format_ax(ax4, label="ACES - 2018", cb=True)
    
    pl.tight_layout()
    pl.savefig(f'{basepath}/papers/continuum_data/figures/SgrB2_2013_comparison.png', 
               bbox_inches='tight', dpi=200)


def main():
    """Main function to generate all figures"""
    # Load data
    sgrb2fh, acesmosaicfh, acesmosaicfeatherfh = load_data()
    
    # Reproject data
    aces_to_sgr, aces_to_sgr_feather = reproject_data(sgrb2fh, acesmosaicfh, acesmosaicfeatherfh)
    
    # Convolve data
    sgr_conv_ACES, SgrB2Beam, ACESbeam, ww = convolve_data(sgrb2fh, acesmosaicfh)
    
    # Calculate conversion factors
    jtok_ACES, jtok_SgrB2 = calculate_jtok(ACESbeam, SgrB2Beam)
    
    # Create the zoomed comparison plot
    print("running plot_comparison_zoom")
    plot_comparison_zoom(sgrb2fh, sgr_conv_ACES, aces_to_sgr, ww, jtok_SgrB2, jtok_ACES, mymap)
    
    # Create the zoomed comparison plot with feathered data
    print("running plot_comparison_zoom w/feather")
    plot_comparison_zoom(sgrb2fh, sgr_conv_ACES, aces_to_sgr_feather, ww, jtok_SgrB2, jtok_ACES, mymap, feathered=True)
    
    # Create the full comparison plot
    print("running plot_comparison_full")
    plot_comparison_full(sgrb2fh, sgr_conv_ACES, aces_to_sgr, ww, jtok_SgrB2, jtok_ACES, mymap)


if __name__ == "__main__":
    main()
