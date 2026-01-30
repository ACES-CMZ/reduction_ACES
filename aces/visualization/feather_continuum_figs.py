import os
from astropy.wcs import WCS
from astropy import wcs
from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.io import fits
from astropy import units as u
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
import matplotlib.colors as mcolors
import numpy as np
import mpl_plot_templates
import mpl_plot_templates.inset_plots
from radio_beam import Beam
import radio_beam
import regions
from astropy import coordinates
from astropy import visualization
import matplotlib.pyplot as pl
import PIL
import pyavm

from aces.visualization.figure_configuration import mymap, mymap2, distance
from aces import conf


# Configuration
basepath = conf.basepath
TENS_DIR = f'{basepath}/TENS/'


# Define regions
b1reg = regions.RectangleSkyRegion(
    center=SkyCoord(0.5*u.deg, -0.06*u.deg, frame='galactic'),
    width=6*u.arcmin,
    height=6*u.arcmin
)


def load_mustang_data():
    """Load MUSTANG data and create necessary objects"""
    fn = 'SgrB2_flux_cut_dt6_filtp05to49_noShift_final_map.fits'
    #mustang_fn = f'{TENS_DIR}/{fn}'
    mustang_plus_planck_fn = f"{TENS_DIR}/{fn.replace('.fits', '')}_PlanckCombined.fits"

    mustang_planck_image = fits.open(mustang_plus_planck_fn)
    mustangdata = mustang_planck_image[0].data
    mustangwcs = wcs.WCS(mustang_planck_image[0].header)
    mustangdata = convolve_fft(mustangdata, Gaussian2DKernel(1), allow_huge=True)
    mustangbeam = radio_beam.Beam.from_fits_header(mustang_planck_image[0].header)

    return mustangdata, mustangwcs, mustangbeam


def load_aces_data():
    """Load ACES data and create necessary objects"""
    ACES = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits')
    ACESheader = ACES[0].header
    ACESwcs = WCS(ACESheader)
    ACESheader['BUNIT'] = 'K'
    ACESbeam = Beam.from_fits_header(ACESheader)
    ACESdata = ACES[0].data
    # 97.3 GHz is roughly the composite effective frequency for alpha=0: see continuumselection / nueff
    jtok_aces = (1*u.Jy).to(u.K, ACESbeam.jtok_equiv(97.3*u.GHz)).value
    ACESdata *= jtok_aces
    return ACESdata, ACESwcs, ACESbeam


def cutout(region, data, wcs):
    """
    Cut out a region from data using the given WCS
    """
    mask = region.to_pixel(wcs).to_mask()
    slices_large, slices_small = mask.get_overlap_slices(data.shape)
    return data[slices_large], wcs[slices_large]


def compare_feathered_to_unfeathered():
    #blcb1 = coordinates.SkyCoord(0.6*u.deg, -0.12*u.deg, frame='galactic')
    #trcb1 = coordinates.SkyCoord(0.44*u.deg, 0.03*u.deg, frame='galactic')

    # Load data
    mustangdata, mustangwcs, mustangbeam = load_mustang_data()
    ACESdata, ACESwcs, ACESbeam = load_aces_data()

    # Load feathered data
    mustangfeathered = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits')
    feathwcs = WCS(mustangfeathered[0].header)

    pl.figure(figsize=(8.5, 3.7*8.5/12), dpi=250)

    # MUSTANG plot
    cmust, wcmust = cutout(b1reg, mustangdata, mustangwcs)
    ax = pl.subplot(1, 3, 1, projection=wcmust)
    im = ax.imshow(cmust, norm=visualization.simple_norm(cmust, stretch='asinh', min_percent=0.5, max_percent=99.95), cmap=mymap)

    ax.set_xlabel("Galactic Longitude")
    ax.set_ylabel("Galactic Latitude")
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.text(0.5, 0.9, "GBT MUSTANG", transform=ax.transAxes, ha='center', color='w')
    cb = pl.colorbar(im)
    cb.set_ticks([0, 0.02, 0.05, 0.1, 0.15, 0.25])

    xell, yell = 10, 10
    ell = mustangbeam.ellipse_to_plot(xell, yell, wcmust.proj_plane_pixel_area()**0.5)
    ell.set_facecolor('none')
    ell.set_edgecolor('red')
    ax.add_artist(ell)

    # ALMA 12m plot
    c12m, wc12m = cutout(b1reg, ACESdata, ACESwcs)
    ax = pl.subplot(1, 3, 3, projection=wc12m)
    im = ax.imshow(c12m,
                   norm=visualization.simple_norm(c12m, stretch='asinh', min_percent=1, max_percent=99.95), cmap='gray_r')
    ax.set_xlabel("Galactic Longitude")
    ax.set_ylabel("Galactic Latitude")
    ax.coords[1].set_ticklabel_visible(False)
    ax.coords[1].set_axislabel("")
    ax.coords[0].set_major_formatter('d.dd')
    ax.text(0.5, 0.9, "ALMA 12m", transform=ax.transAxes, ha='center')
    pl.colorbar(im)

    xell, yell = 30, 30
    ell = ACESbeam.ellipse_to_plot(xell, yell, wc12m.proj_plane_pixel_area()**0.5)
    ell.set_facecolor('none')
    ell.set_edgecolor('red')
    ax.add_artist(ell)

    # Feathered plot
    fma, wfma = cutout(b1reg, mustangfeathered[0].data, feathwcs)
    ax = pl.subplot(1, 3, 2, projection=wfma)
    im = ax.imshow(fma,
                   norm=visualization.simple_norm(fma, stretch='asinh', min_percent=0.5, max_percent=99.95), cmap=mymap)
    ax.set_xlabel("Galactic Longitude")
    ax.set_ylabel("Galactic Latitude")
    ax.coords[1].set_ticklabel_visible(False)
    ax.coords[1].set_axislabel("")
    ax.coords[0].set_major_formatter('d.dd')
    ax.text(0.5, 0.9, "Feathered", transform=ax.transAxes, ha='center', color='w')
    cb = pl.colorbar(im)
    cb.set_ticks([0, 0.02, 0.05, 0.1, 0.15, 0.25])

    ell = ACESbeam.ellipse_to_plot(xell, yell, wc12m.proj_plane_pixel_area()**0.5)
    ell.set_facecolor('none')
    ell.set_edgecolor('red')
    ax.add_artist(ell)

    for ax in pl.gcf().axes:
        try:
            ax.coords[0].set_ticklabel(rotation=45, pad=30)
        except Exception as ex:  # noqa: F841
            pass

    pl.tight_layout()
    pl.subplots_adjust(hspace=0, wspace=0.0)  # adjust by changing figsize
    pl.savefig("/orange/adamginsburg/web/secure/ACES/mosaics/MUSTANG_vs_ACES_zoomSgrB1.pdf", bbox_inches='tight')
    pl.savefig("/orange/adamginsburg/web/secure/ACES/mosaics/MUSTANG_vs_ACES_zoomSgrB1.png", bbox_inches='tight')


def mustang_feather_zoomregions():
    # Load ACES data
    ACESdata, ACESwcs, ACESbeam = load_aces_data()

    # Load feathered data
    mustangfeatheredfn = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'
    rslt = fits.open(mustangfeatheredfn)[0]
    norm = simple_norm(rslt.data[:, :], vmin=0.0001, vmax=1.5, stretch='log')
    colordata = mymap(norm(rslt.data[:, :]))

    # flip the Y-axis
    ct = (colordata[::-1, :, :3] * 256)
    ct[ct > 255] = 255
    img = PIL.Image.fromarray(ct.astype('uint8'))
    img.save(f'{basepath}/mosaics/continuum/MUSTANG_12m_feather_noaxes.png')

    avm = pyavm.AVM.from_wcs(ACESwcs)
    avm.embed(f'{basepath}/mosaics/continuum/MUSTANG_12m_feather_noaxes.png',
              f'{basepath}/mosaics/continuum/MUSTANG_12m_feather_noaxes.png')

    from continuum_fig1a import plot_fullwidth_figure

    fig = plot_fullwidth_figure(mustangfeatheredfn, scale=1, normkwargs={'vmin': 0.0001, 'vmax': 1.5, 'stretch': 'log'},
                                colormap=mymap,
                                label=r"$T_B$ [K]", figsize=(14, 4), dpi=250, scalebar=True, beam=True,
                                scalebar_color='w',
                                subset_label=False)

    cb = fig.axes[0].images[-1].colorbar
    cb.set_ticks([0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.8, 1.4])

    fig.savefig(f'{basepath}/mosaics/continuum/MUSTANG_12m_feather.png', bbox_inches='tight', dpi=300)
    fig.savefig(f'{basepath}/mosaics/continuum/MUSTANG_12m_feather.pdf', bbox_inches='tight', dpi=300)

    # zoomreg overlay
    zoomregs = regions.Regions.read(f'{basepath}/reduction_ACES/aces/data/regions/zoom_regions.reg')

    ax = fig.gca()

    for reg in zoomregs:
        name = reg.meta['text']
        preg = reg.to_pixel(ACESwcs)
        rect = preg.plot(ax=ax)
        rect.set_linewidth(0.5)

        xr = rect.get_x()
        yr = rect.get_y()
        width = rect.get_width()
        height = rect.get_height()
        text_x = xr + width / 2  # Center horizontally
        if name in ("Pistol Sickle Quintuplet", "20 km/s Cloud", "G359.63"):
            text_y = yr - 0.1
            va = 'top'
        else:
            text_y = yr + height + 0.1  # Slightly above the rectangle
            va = 'bottom'
        ax.text(text_x, text_y, name, ha='center', va=va, color=rect.get_edgecolor(), fontsize=8)

        ax.coords[0].set_major_formatter('d.dd')
        ax.coords[1].set_major_formatter('d.dd')
        ax.coords[0].set_ticks(spacing=0.25*u.deg)

    fig.savefig(f'{basepath}/mosaics/continuum/MUSTANG_12m_feather_zoomregions.png', bbox_inches='tight', dpi=300)

    # afmhot
    fig = plot_fullwidth_figure(mustangfeatheredfn, scale=1, normkwargs={'vmin': 0.0001, 'vmax': 1.5, 'stretch': 'log'},
                                colormap=mymap2,
                                label=r"$T_B$ [K]", figsize=(14, 4), dpi=250, scalebar=True, beam=True,
                                scalebar_color='w',
                                subset_label=False)

    cb = fig.axes[0].images[-1].colorbar
    cb.set_ticks([0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.8, 1.4])

    fig.savefig(f'{basepath}/mosaics/continuum/MUSTANG_12m_feather_afmhot.png', bbox_inches='tight', dpi=300)

    TENS_DIR = '/orange/adamginsburg/ACES/TENS/'
    mustang_name = 'SgrB2_flux_cut_dt6_filtp05to49_noShift_final_map_PlanckCombined.fits'
    mustangfn = f'{TENS_DIR}/{mustang_name}'

    fig = plot_fullwidth_figure(mustangfn, scale=1, normkwargs={'vmin': 0.0001, 'vmax': 1.5, 'stretch': 'log'},
                                colormap=mymap,
                                label=r"$T_B$ [K]", figsize=(14, 4), dpi=250, scalebar=True, beam=True,
                                scalebar_color='w',
                                scalebarstart=SkyCoord(-0.75*u.deg, -0.4*u.deg, frame='galactic'),
                                subset_label=False)

    cb = fig.axes[0].images[-1].colorbar
    cb.set_ticks([0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.8, 1.4])
    
    fig.axes[0].coords[0].set_ticks(spacing=0.5*u.deg)

    fig.savefig(f'{basepath}/mosaics/continuum/MUSTANG.png', bbox_inches='tight', dpi=300)
    fig.savefig(f'{basepath}/mosaics/continuum/MUSTANG.pdf', bbox_inches='tight', dpi=300)



if __name__ == '__main__':
    mustang_feather_zoomregions()
    compare_feathered_to_unfeathered()