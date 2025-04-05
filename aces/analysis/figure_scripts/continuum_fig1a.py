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

from aces.analysis.figure_scripts.figure_configuration import mymap, mymap2, pl, distance

from aces import conf
basepath = conf.basepath
diag_dir = f'{basepath}/diagnostic_plots/'


for colormap, cmname in ((mymap, 'grey-hot'), (mymap2, 'grey-afmhot')):
    
    for fn, name in zip((f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits',
                         f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits',
                         f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_mosaic.fits',
                        ),
                        ('aggregate', 'high', 'low')) :
        fig = pl.figure(figsize=(14,4), dpi=150)
        imgdata = fits.open(fn)[0].data
        imgdata[imgdata == 0] = np.nan
        rmswcs = WCS(fits.open(fn)[0].header)
        ax = pl.subplot(projection=rmswcs)
        #im = ax.imshow(imgdata, norm=simple_norm(imgdata, min_cut=-0.0001, max_cut=0.0006, stretch='linear'),
        #          cmap='gray_r')
        #im2 = ax.imshow(imgdata, norm=simple_norm(imgdata, min_cut=0.0006, max_percent=99.9, stretch='log'),
        #          cmap='inferno')
        im = ax.imshow(imgdata*1e3, norm=simple_norm(imgdata*1e3, vmin=-0.5, vmax=20, stretch='asinh'),
                  cmap=colormap)
        cb = pl.colorbar(mappable=im, pad=0.01);
        cb.set_label(r"$S_\nu$ [mJy/beam]")
        ax.set_xlabel("Galactic Longitude")
        ax.set_ylabel("Galactic Latitude")
        
        ax.coords[0].set_major_formatter('d.dd')
        ax.coords[1].set_major_formatter('d.dd')
        ax.coords[0].set_ticks(spacing=0.25*u.deg)
    
        # ax.text(0.99, 0.93, {'aggregate': 'Aggregate',
        #                      'low': 'spw 25+27',
        #                      'high': 'spw 33+35'}[name],
        #         transform=ax.transAxes,
        #         horizontalalignment='right')
        
        ax.coords[0].set_major_formatter('d.dd')
        ax.coords[1].set_major_formatter('d.dd')
        ax.coords[0].set_ticks(spacing=0.25*u.deg)
        
        mpl_plot_templates.inset_plots.make_scalebar(ax, left_side=SkyCoord(0.5*u.deg, -0.28*u.deg, frame='galactic'),
                                                     length=(25*u.pc/distance).to(u.arcsec, u.dimensionless_angles()),
                                                     color='k', label='25 pc')
        
        
        
        fig.savefig(f'{diag_dir}/FullField12m_circularbeam_{name}_{cmname}.png', bbox_inches='tight', dpi=250)