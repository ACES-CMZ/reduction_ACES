from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
import matplotlib.colors as mcolors
import numpy as np
import mpl_plot_templates
import mpl_plot_templates.inset_plots
from mpl_toolkits.axes_grid1 import make_axes_locatable

from aces.visualization.figure_configuration import mymap, mymap2, pl, distance
from aces.visualization.compare_to_other import format_ax

import radio_beam

from aces import conf
basepath = conf.basepath
diag_dir = f'{basepath}/diagnostic_plots/'


def plot_fullwidth_figure(fn, imgdata=None,
                          name='Aggregate',
                          label=r"$S_\nu$ [mJy/beam]",
                          figsize=(14, 4), dpi=250, scalebar=True, beam=True,
                          colormap=mymap,
                          normkwargs={'vmin': -0.5, 'vmax': 20, 'stretch': 'asinh'},
                          scale=1e3,
                          scalebar_color='k',
                          subset_label=False,
                          ax=None,
                          fig=None):
    if fig is None:
        fig = pl.figure(figsize=figsize, dpi=dpi)
    if imgdata is None:
        imgdata = fits.open(fn)[0].data
        imgdata[imgdata == 0] = np.nan
    ww = WCS(fits.open(fn)[0].header)
    if ax is None:
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

def plot_three_panel_figure(fn1, fn2, fn3, names=None, colormap=None,
                            label=r"$S_\nu$ [mJy/beam]",
                            dpi=300,
                            subset_label=True, scale=1e3, normkwargs=None,
                            scalebar=True,
                            scalebar_color='k',
                            cbar_size="2%"):
    """
    Create a 3-panel figure with vertically stacked panels sharing a common X-axis.
    
    Parameters
    ----------
    fn1, fn2, fn3 : str
        Filenames for the three panels
    names : list of str, optional
        Names for each panel
    colormap : str or colormap object
        Colormap to use for all panels
    label : str
        Label for the colorbar
    subset_label : bool
        Whether to add subset labels
    scale : float
        Scaling factor for the data
    normkwargs : dict
        Keyword arguments for the normalization
    """

    with pl.rc_context({'font.size': 12}):
        # Create figure with no spacing between subplots
        fig = pl.figure(figsize=(8, 7), dpi=dpi)
        fig.subplots_adjust(hspace=0, wspace=0)
        
        # Create three subplots sharing the x-axis with GridSpec for tight control
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(3, 1, height_ratios=[1, 1, 1], hspace=0)
        
        # Create three subplots sharing the x-axis
        ax1 = fig.add_subplot(gs[0], projection=WCS(fits.open(fn1)[0].header))
        ax2 = fig.add_subplot(gs[1], sharex=ax1, projection=WCS(fits.open(fn2)[0].header))
        ax3 = fig.add_subplot(gs[2], sharex=ax1, projection=WCS(fits.open(fn3)[0].header))
        
        # Use default normkwargs if none provided
        if normkwargs is None:
            normkwargs = {'vmin': -0.5, 'vmax': 20, 'stretch': 'asinh'}
        
        # Plot each panel with matched colorbar heights
        # Handle the top panel (ax1)
        imgdata1 = fits.open(fn1)[0].data
        if names and names[0]:
            name1 = names[0]
        else:
            name1 = None
        imgdata1[imgdata1 == 0] = np.nan
        im1 = ax1.imshow(imgdata1*scale, norm=simple_norm(imgdata1*scale, **normkwargs),
                    cmap=colormap)
        cb1 = format_ax(ax1, label=label, hidex=True, cbar_size=cbar_size)
        ax1.set_ylabel("Galactic Latitude")
        ax1.tick_params(axis='x', labelbottom=False)  # Hide x labels for top panel
        
        # Handle the middle panel (ax2)
        imgdata2 = fits.open(fn2)[0].data
        if names and len(names) > 1 and names[1]:
            name2 = names[1]
        else:
            name2 = None
        imgdata2[imgdata2 == 0] = np.nan
        im2 = ax2.imshow(imgdata2*scale, norm=simple_norm(imgdata2*scale, **normkwargs),
                    cmap=colormap)
        cb2 = format_ax(ax2, label=label, hidex=True, cbar_size=cbar_size)
        ax2.set_ylabel("Galactic Latitude")
        ax2.tick_params(axis='x', labelbottom=False)  # Hide x labels for middle panel
        
        # Handle the bottom panel (ax3)
        imgdata3 = fits.open(fn3)[0].data
        if names and len(names) > 2 and names[2]:
            name3 = names[2]
        else:
            name3 = None
        imgdata3[imgdata3 == 0] = np.nan
        im3 = ax3.imshow(imgdata3*scale, norm=simple_norm(imgdata3*scale, **normkwargs),
                    cmap=colormap)
        cb3 = format_ax(ax3, label=label, cbar_size=cbar_size)
        ax3.set_xlabel("Galactic Longitude")
        ax3.set_ylabel("Galactic Latitude")

        labeldict = {
            'aggregate': 'Aggregate',
            'low': 'spw 25+27',
            'high': 'spw 33+35',
            'spw25_27': 'spw 25+27',
            'spw33_35': 'spw 33+35',
        }
        
        # Add subset labels if requested
        if subset_label and names:
            if name1:
                ax1.text(0.99, 0.93, labeldict.get(name1, name1),
                        transform=ax1.transAxes,
                        horizontalalignment='right')
            if name2:
                ax2.text(0.99, 0.93, labeldict.get(name2, name2),
                        transform=ax2.transAxes,
                        horizontalalignment='right')
            if name3:
                ax3.text(0.99, 0.93, labeldict.get(name3, name3),
                        transform=ax3.transAxes,
                        horizontalalignment='right')
                        
        # Configure axis formatting
        for ax in [ax1, ax2, ax3]:
            ax.coords[0].set_major_formatter('d.dd')
            ax.coords[1].set_major_formatter('d.dd')
            ax.coords[0].set_ticks(spacing=0.25*u.deg)

        if scalebar:
            for ax in [ax1, ax2, ax3]:
                mpl_plot_templates.inset_plots.make_scalebar(ax,
                                                            left_side=SkyCoord(0.5*u.deg, -0.28*u.deg, frame='galactic'),
                                                            length=(25*u.pc/distance).to(u.arcsec, u.dimensionless_angles()),
                                                            color=scalebar_color,
                                                            label='25 pc')

        
        # Remove the space between the plots
        plt = pl
        plt.subplots_adjust(hspace=0)
        
        # Make sure no white space is visible between subplots
        for ax in [ax1, ax2]:
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.spines['bottom'].set_visible(False)
            
        for ax in [ax2, ax3]:
            ax.spines['top'].set_visible(False)
            
        # Manually adjust subplot positions to completely eliminate vertical space
        # This gets the current positions of the subplots
        pos1 = ax1.get_position()
        pos2 = ax2.get_position()
        pos3 = ax3.get_position()
        
        # Calculate new positions to make them touch without any gap
        # Keep the same width and x position
        height = pos1.height  # Assuming all heights are equal
        
        # Set new positions - direct bottom-to-top alignment
        ax1.set_position([pos1.x0, pos3.y0 + 2*height, pos1.width, height])
        ax2.set_position([pos2.x0, pos3.y0 + height, pos2.width, height])
        # ax3 stays at its position
        
        # Execute tight_layout again to adjust margins properly
        # But exclude the colorbar axes from consideration
        colorbar_axes = [cb1.ax, cb2.ax, cb3.ax]
        all_axes = fig.get_axes()
        for ax in colorbar_axes:
            if ax in all_axes:
                all_axes.remove(ax)
        
        # Apply tight_layout to the remaining axes
        fig.tight_layout(h_pad=0, rect=[0, 0, 0.95, 1])
        
        return fig


if __name__ == '__main__':
    for colormap, cmname in ((mymap, 'grey-hot'), (mymap2, 'grey-afmhot')):

        filenames =(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits',
                    f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits',
                    f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_mosaic.fits',)
        fname_labels = ('aggregate', 'high', 'low')
        iterables = zip(filenames, fname_labels)
        for fn, name in iterables:
            fig = plot_fullwidth_figure(fn, name=name, colormap=colormap)
            fig.savefig(f'{diag_dir}/FullField12m_circularbeam_{name}_{cmname}.png', bbox_inches='tight', dpi=250)
            fig.savefig(f'{diag_dir}/FullField12m_circularbeam_{name}_{cmname}.pdf', bbox_inches='tight', dpi=250)

        fig = plot_three_panel_figure(*filenames, names=fname_labels, colormap=colormap, subset_label=True)
        fig.savefig(f'{diag_dir}/FullField12m_circularbeam_three_panel_{cmname}.png', bbox_inches='tight', dpi=300)
        fig.savefig(f'{diag_dir}/FullField12m_circularbeam_three_panel_{cmname}.pdf', bbox_inches='tight', dpi=300)

        iterables = (
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_maskedrms_mosaic.fits', 'aggregate'),
            #(f'{basepath}/mosaics/continuum/12m_continuum_residual_commonbeam_circular_reimaged_maskedrms_mosaic.fits', 'fullresidual'),
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_maskedrms_mosaic.fits', 'spw33_35'),
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_maskedrms_mosaic.fits', 'spw25_27'),
        )
        for fn, name in iterables:
            fig = plot_fullwidth_figure(fn, name=name, label="RMS [mJy beam$^{-1}$]", subset_label=True, colormap=colormap)
            fig.savefig(f'{diag_dir}/RMS_map_{name}_{cmname}.png', bbox_inches='tight', dpi=250)

        filenames = [x[0] for x in iterables]
        fname_labels = [x[1] for x in iterables]
        fig = plot_three_panel_figure(*filenames, names=fname_labels, colormap=colormap, subset_label=True,
                                      normkwargs = {'vmin': 0.0, 'vmax': 1, 'stretch': 'asinh'})
        fig.savefig(f'{diag_dir}/RMS_three_panel_{cmname}.png', bbox_inches='tight', dpi=300)
        fig.savefig(f'{diag_dir}/RMS_three_panel_{cmname}.pdf', bbox_inches='tight', dpi=300)


        iterables = (
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_maskedrms_mosaic.fits',
             f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits', 'aggregate'),
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_maskedrms_mosaic.fits',
             f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits', 'spw33_35'),
            (f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_maskedrms_mosaic.fits',
             f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_mosaic.fits', 'spw25_27'),
        )

        snmaps = []
        for rmsfn, fn, name in iterables:

            map = fits.open(fn)[0].data
            rmsmap = fits.open(rmsfn)[0].data
            snmap = map / rmsmap
            snmap[rmsmap == 0] = np.nan
            snmap[snmap == 0] = np.nan
            snmaps.append((fn, snmap, name))
            fig = plot_fullwidth_figure(fn, imgdata=snmap, name=name, label="S/N Ratio", subset_label=True, colormap=colormap, scale=1)
            fig.savefig(f'{diag_dir}/SignalToNoise_map_{name}_{cmname}.png', bbox_inches='tight', dpi=250)
            fig.savefig(f'{diag_dir}/SignalToNoise_map_{name}_{cmname}.pdf', bbox_inches='tight', dpi=250)
        
        
        # Create three panel figure for S/N maps
        fig = plot_three_panel_figure(
            snmaps[0][0], snmaps[1][0], snmaps[2][0],
            names=[snmaps[0][2], snmaps[1][2], snmaps[2][2]],
            colormap=colormap, 
            label="S/N Ratio",
            subset_label=True,
            scale=1,
            normkwargs={'vmin': 0, 'vmax': 50, 'stretch': 'linear'}
        )
        fig.savefig(f'{diag_dir}/SignalToNoise_three_panel_{cmname}.png', bbox_inches='tight', dpi=300)
        fig.savefig(f'{diag_dir}/SignalToNoise_three_panel_{cmname}.pdf', bbox_inches='tight', dpi=300)
