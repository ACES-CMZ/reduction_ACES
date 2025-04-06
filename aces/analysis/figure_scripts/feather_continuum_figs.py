from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.visualization import simple_norm
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
import matplotlib.colors as mcolors
import numpy as np
import mpl_plot_templates
import mpl_plot_templates.inset_plots
from radio_beam import Beam
import regions

from aces.analysis.figure_scripts.figure_configuration import mymap, mymap2, pl, distance, distance

import PIL
import pyavm

from aces import conf
basepath = conf.basepath


if __name__ == '__main__':

    ACES = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits')
    ACESheader = ACES[0].header
    ACESwcs = WCS(ACES[0].header)
    ACESheader['BUNIT'] = 'K'
    ACESbeam = Beam.from_fits_header(ACESheader)
    ACESdata = ACES[0].data
    # 97.3 GHz is roughly the composite effective frequency for alpha=0: see continuumselection / nueff
    jtok_aces = (1*u.Jy).to(u.K, Beam.from_fits_header(ACESheader).jtok_equiv(97.3*u.GHz)).value
    # DEBUG: added *10
    ACESdata *= jtok_aces

    mustangfeatheredfn = f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'
    rslt = fits.open(mustangfeatheredfn)[0]
    norm = simple_norm(rslt.data[:,:], vmin=0.0001, vmax=1.5, stretch='log')
    colordata = mymap(norm(rslt.data[:,:]))

    # flip the Y-axis
    ct = (colordata[::-1,:,:3] * 256)
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

