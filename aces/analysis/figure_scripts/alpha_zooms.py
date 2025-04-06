from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
from spectral_cube import SpectralCube
import spectral_cube
import matplotlib.colors as mcolors
import numpy as np
import mpl_plot_templates
from radio_beam import Beam
import regions
from scipy import ndimage

from aces.analysis.figure_scripts.figure_configuration import mymap, mymap2, mymap_gray, pl, distance

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib import pyplot as plt

from aces import conf
basepath = conf.basepath


def add_cb(ax, label):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cax.set_xticks([])
    cax.set_xlabel("")
    cax.set_xticklabels([])
    im = ax.images[-1]
    cb = plt.colorbar(im, cax=cax)
    
    cax.coords[0].set_axislabel('')
    cax.coords[0].set_ticklabel_visible(False)
    cax.coords[0].set_ticks_visible(False)
    cax.coords[1].set_ticks_position('r')
    cax.coords[1].set_ticklabel_position('r')
    cax.coords[1].set_axislabel_position('r')
    cax.coords[1].set_axislabel(label)
    return cb


pl.rcParams['font.size'] = 14

zoomregs = regions.Regions.read(f'{basepath}/reduction_ACES/aces/data/regions/zoom_regions.reg')

manual_alpha = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_manual_alpha_mosaic.fits')
tt_alpha = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_alpha_mosaic.fits')
#tt_alpha_err = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_alphaerror_mosaic.fits')

lofrqhdu = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_mosaic.fits')
hifrqhdu = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits')
lofrq = spectral_cube.Projection.from_hdu(lofrqhdu)
hifrq = spectral_cube.Projection.from_hdu(hifrqhdu)

lofrqrmshdu = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_maskedrms_mosaic.fits')
hifrqrmshdu = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_maskedrms_mosaic.fits')
lofrqrms = spectral_cube.Projection.from_hdu(lofrqrmshdu)
hifrqrms = spectral_cube.Projection.from_hdu(hifrqrmshdu)

snr_lo = lofrq / lofrqrms

featherfh = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits')
feather = spectral_cube.Projection.from_hdu(featherfh)

lothresh = snr_lo > 0.5
hithresh = snr_lo > 6
joint_mask = ndimage.binary_dilation(hithresh, iterations=-1, mask=lothresh)

for reg in zoomregs:
    preg = reg.to_pixel(lofrq.wcs)
    mask = preg.to_mask()
    slc_lg, slc_sm = mask.get_overlap_slices(lofrq.shape)
    ww_sl = lofrq.wcs[slc_lg]

    lo_cutout = lofrq[slc_lg]
    xell, yell = np.array(lo_cutout.shape) * 0.05

    pl.figure(figsize=(25, 4))
    ax1 = pl.subplot(1,4,1, projection=ww_sl)
    norm = simple_norm(lo_cutout.value*1e3, min_percent=0.5, max_percent=99.95)
    im = ax1.imshow(lo_cutout.value*1e3, cmap=mymap, norm=norm)
    ell = lofrq.beam.ellipse_to_plot(xell, yell, ww_sl.proj_plane_pixel_area()**0.5)
    ell.set_facecolor('none')
    ell.set_edgecolor('red')
    ax1.add_artist(ell)                               
    cb = add_cb(ax1, "S$_{25+27}$ [mJy beam$^{-1}$]")
    
    
    ax2 = pl.subplot(1,4,2, projection=ww_sl)
    im = ax2.imshow(hifrq[slc_lg].value*1e3, cmap=mymap, norm=norm)
    ell = hifrq.beam.ellipse_to_plot(xell, yell, ww_sl.proj_plane_pixel_area()**0.5)
    ell.set_facecolor('none')
    ell.set_edgecolor('red')
    ax2.add_artist(ell)                               
    cb = add_cb(ax2, "S$_{33+35}$ [mJy beam$^{-1}$]")
    
    ax3 = pl.subplot(1,4,3, projection=ww_sl)
    alpha_cutout = manual_alpha[0].data[slc_lg]
    bad_mask = np.ones_like(alpha_cutout)
    bad_mask[~joint_mask[slc_lg]] = np.nan
    im0 = ax3.imshow(hifrq[slc_lg].value*1e3, cmap=mymap_gray, norm=norm, zorder=-5)
    im3 = ax3.imshow(alpha_cutout * bad_mask, vmin=-5, vmax=5, cmap='Spectral')
    cm = add_cb(ax3, r"Spectral Index $\alpha$",)
    
    ax4 = pl.subplot(1,4,4, projection=ww_sl)
    im = ax4.imshow(feather.value[slc_lg], cmap=mymap,
                    norm=simple_norm(feather.value[slc_lg], min_percent=0.5, max_percent=99.95))
    cm = add_cb(ax4, r"T$_B$ [K]")

    name = reg.meta['text']

    for ax in (ax2, ax3, ax4):
        lat = ax.coords[1]
        lat.set_axislabel('')
        lat.set_ticklabel_visible(False)
    ax1.coords[1].set_axislabel("Galactic Latitude")
    for ax in (ax1, ax2, ax3, ax4):
        lon, lat = ax.coords
        lon.set_axislabel("Galactic Longitude")
        lon.set_major_formatter('d.dd')
        lat.set_major_formatter('d.dd')
        if name in ('Sgr C', '20 km/s Cloud', 'Sgr A*'):
            lon.set_ticklabel(rotation=45, pad=40)
    pl.suptitle(name)
    if name in ('North Arches', '20 km/s Cloud'):
        pl.subplots_adjust(wspace=0.2, hspace=0.01)
    elif name == 'Sgr A North':
        pl.subplots_adjust(wspace=-0.1, hspace=0.01)
    else:
        pl.subplots_adjust(wspace=0.03, hspace=0.01)

    pl.savefig(f'{basepath}/papers/continuum_data/figures/{name.replace(" ", "_").replace("/","").replace("*","star")}_hi_lo_alpha.pdf', bbox_inches='tight')