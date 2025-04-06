from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
import matplotlib.colors as mcolors
import numpy as np
import mpl_plot_templates
import mpl_plot_templates.inset_plots

from aces.analysis.figure_scripts.figure_configuration import mymap, mymap2, pl, distance

import radio_beam

from aces import conf
basepath = conf.basepath
diag_dir = f'{basepath}/diagnostic_plots/'


def plot_fullwidth_figure(fn, imgdata=None,
                          name='Aggregate',
                          label=r"$S_\nu$ [mJy/beam]", figsize=(14, 4), dpi=250, scalebar=True, beam=True,
                          colormap=mymap,
                          normkwargs={'vmin': -0.5, 'vmax': 20, 'stretch': 'asinh'},
                          scale=1e3,
                          scalebar_color='k',
                          subset_label=False):
    fig = pl.figure(figsize=figsize, dpi=dpi)
    if imgdata is None:
        imgdata = fits.open(fn)[0].data
        imgdata[imgdata == 0] = np.nan
    ww = WCS(fits.open(fn)[0].header)
    ax = pl.subplot(projection=ww)
    #im = ax.imshow(imgdata, norm=simple_norm(imgdata, min_cut=-0.0001, max_cut=0.0006, stretch='linear'),
    #          cmap='gray_r')
    #im2 = ax.imshow(imgdata, norm=simple_norm(imgdata, min_cut=0.0006, max_percent=99.9, stretch='log'),
    #          cmap='inferno')
    im = ax.imshow(imgdata*scale, norm=simple_norm(imgdata*scale, **normkwargs),
                   cmap=colormap)
    cb = pl.colorbar(mappable=im, pad=0.01)
    cb.set_label(label)
    ax.set_xlabel("Galactic Longitude")
    ax.set_ylabel("Galactic Latitude")

    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.coords[0].set_ticks(spacing=0.25*u.deg)

    if subset_label:
        ax.text(0.99, 0.93, {'aggregate': 'Aggregate',
                             'low': 'spw 25+27',
                             'high': 'spw 33+35',
                             'spw25_27': 'spw 25+27',
                             'spw33_35': 'spw 33+35',
                             }[name],
                transform=ax.transAxes,
                horizontalalignment='right')

    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.coords[0].set_ticks(spacing=0.25*u.deg)

    if scalebar:
        mpl_plot_templates.inset_plots.make_scalebar(ax,
                                                     left_side=SkyCoord(0.5*u.deg, -0.28*u.deg, frame='galactic'),
                                                     length=(25*u.pc/distance).to(u.arcsec, u.dimensionless_angles()),
                                                     color=scalebar_color,
                                                     label='25 pc')

    if beam:
        beam = radio_beam.Beam.from_fits_header(fn)
        pixscale = (ww.proj_plane_pixel_area()**0.5).to(u.arcsec)
        ellipse_artist = beam.ellipse_to_plot(50, 50, pixscale)
        ellipse_artist.set_facecolor('b')
        ellipse_artist.set_edgecolor('b')
        ellipse_artist.set_linewidth(0.25)
        ellipse_artist.set_fill(True)
        ax.add_artist(ellipse_artist)
        rectangle_artist = pl.Rectangle((20, 20,), 60, 60, edgecolor='k', facecolor='none', lw=0.25)
        ax.add_artist(rectangle_artist)

    return fig


if __name__ == '__main__':
    for colormap, cmname in ((mymap, 'grey-hot'), (mymap2, 'grey-afmhot')):

        for fn, name in zip((f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits',
                             f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits',
                             f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_mosaic.fits',
                            ), ('aggregate', 'high', 'low')):
            fig = plot_fullwidth_figure(fn, name=name, colormap=colormap)
            fig.savefig(f'{diag_dir}/FullField12m_circularbeam_{name}_{cmname}.png', bbox_inches='tight', dpi=250)

        iterables = (
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_maskedrms_mosaic.fits', 'aggregate'),
            #(f'{basepath}/mosaics/continuum/12m_continuum_residual_commonbeam_circular_reimaged_maskedrms_mosaic.fits', 'fullresidual'),
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_maskedrms_mosaic.fits', 'spw33_35'),
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_maskedrms_mosaic.fits', 'spw25_27'),
        )
        for fn, name in iterables:
            fig = plot_fullwidth_figure(fn, name=name, label="RMS [mJy beam$^{-1}$]", subset_label=True, colormap=colormap)
            fig.savefig(f'{diag_dir}/RMS_map_{name}_{cmname}.png', bbox_inches='tight', dpi=250)

        iterables = (
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_maskedrms_mosaic.fits',
             f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits', 'aggregate'),
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_maskedrms_mosaic.fits',
             f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits', 'spw33_35'),
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_maskedrms_mosaic.fits',
             f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_mosaic.fits', 'spw25_27'),
        )

        for rmsfn, fn, name in iterables:

            map = fits.open(fn)[0].data
            rmsmap = fits.open(rmsfn)[0].data
            snmap = map / rmsmap
            snmap[rmsmap == 0] = np.nan
            snmap[snmap == 0] = np.nan
            fig = plot_fullwidth_figure(fn, imgdata=snmap, name=name, label="S/N Ratio", subset_label=True, colormap=colormap)
            fig.savefig(f'{diag_dir}/SignalToNoise_map_{name}_{cmname}.png', bbox_inches='tight')