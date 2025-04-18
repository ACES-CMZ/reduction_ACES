#!/usr/bin/env python
# coding: utf-8

import numpy as np
import aces
from astropy.table import Table
import pylab as pl
import matplotlib.patches as mpatches
from astropy.io import fits
import radio_beam
from spectral_cube import SpectralCube
from astropy import units as u
from astropy.convolution import Gaussian2DKernel
from astropy.wcs import WCS
from aces.imaging.make_mosaic import rms_map
from aces.visualization.figure_configuration import mymap, mymap2, mymap_gray, pl, distance
from spectral_cube.lower_dimensional_structures import Projection, Slice
import copy
import warnings
from aces import conf
basepath = conf.basepath

diag_dir = f'{basepath}/diagnostic_plots/'


def load_and_filter_data(filename=f'{basepath}/tables/metadata_image.tt0.ecsv'):
    """
    Load data from the specified ECSV file and filter out undesired rows.
    """
    data = Table.read(filename)
    good = ((np.char.find(data['filename'], 'preQA3') == -1) |
            (np.char.find(data['filename'], 'X15a0_X178') == -1) |
            (np.char.find(data['filename'], 'obsolete') == -1))
    data = data[good]
    assert len(data) > 0
    return data


def plot_beam_size_histograms(data):
    """
    Generate histograms of beam sizes for high and aggregate frequencies.
    """
    TM = data['array'] == '12M'
    bins = np.linspace(1.0, 2.2, 10)

    high = np.array([('spw33_35' in fn) and ('v1' not in fn) for fn in data['filename']])
    agg = np.array([('spw25_27_29_31_33_35' in fn) and ('v1' not in fn) for fn in data['filename']])

    pl.figure()
    pl.hist(data['bmaj'][TM & high], bins=bins, label='$B_{maj}$', histtype='step', linewidth=2)
    pl.hist(data['bmin'][TM & high], bins=bins, label='$B_{min}$', histtype='step', linewidth=2, alpha=0.5)
    pl.hist((data['bmin'] * data['bmaj'])[TM & high]**0.5, bins=bins, label='$B_{avg}$', alpha=0.5,
            facecolor='none', edgecolor='k', histtype='step', linewidth=4)
    pl.legend()
    pl.xlabel("Beam FWHM [arcsec]")
    pl.ylabel("Number of Fields")
    pl.savefig(f'{diag_dir}/beam_size_histogram_high.png')

    pl.figure()
    pl.hist(data['bmaj'][TM & agg], bins=bins, label='$B_{maj}$', histtype='step', linewidth=2)
    pl.hist(data['bmin'][TM & agg], bins=bins, label='$B_{min}$', histtype='step', linewidth=2, alpha=0.5)
    pl.hist((data['bmin'] * data['bmaj'])[TM & agg]**0.5, bins=bins, label='$B_{avg}$', alpha=0.5,
            facecolor='none', edgecolor='k', linewidth=4, histtype='step')
    pl.legend()
    pl.xlabel("Beam FWHM [arcsec]")
    pl.ylabel("Number of Fields")
    pl.savefig(f'{diag_dir}/beam_size_histogram_allwindows.png')


def plot_noise_vs_beam(data):
    """
    Generate scatter plot of noise level vs beam geometric average FWHM.
    """
    TM = data['array'] == '12M'
    high = np.array([('spw33_35' in fn) and ('v1' not in fn) for fn in data['filename']])
    low = np.array([('spw25_27.' in fn) and ('v1' not in fn) for fn in data['filename']])
    agg = np.array([('spw25_27_29_31_33_35' in fn) and ('v1' not in fn) for fn in data['filename']])

    pl.figure()
    pl.scatter((data['bmaj'] * data['bmin'])[TM & high]**0.5, data['mad'][TM & high] * 1e3, label='spw33 & 35: 3.6 GHz',
               edgecolor='blue', facecolor='none', linewidth=2)
    pl.scatter((data['bmaj'] * data['bmin'])[TM & low]**0.5, data['mad'][TM & low] * 1e3, label='spw25 & 27: 3.6 GHz',
               marker='s', edgecolor='orange', facecolor='none', linewidth=2)
    pl.scatter((data['bmaj'] * data['bmin'])[TM & agg]**0.5, data['mad'][TM & agg] * 1e3, label='all spw: 4.8 GHz',
               marker="^", facecolor='none', edgecolor='red', linewidth=2)
    pl.axvline(1.5, color='k', linestyle='--')
    pl.text(1.47, 0.25, "1.5\" requested", rotation=90)
    pl.axvline(1.8, color='k', linestyle='--')
    pl.text(1.77, 0.25, "1.8\" ALMA tolerance", rotation=90)
    pl.axhline(0.13, color='k', linestyle='--')
    pl.text(2.13, 0.05, "0.13 mJy requested", horizontalalignment='right', backgroundcolor='w')
    pl.annotate('', (2, 0.13), (1.95, 0.07), arrowprops={'arrowstyle': '->'})
    yl = pl.ylim()
    pl.ylim(0, yl[1])
    pl.legend(fancybox=True, framealpha=0.5)
    pl.xlabel("Beam Geometric Average FWHM [arcsec]")
    pl.ylabel("Noise Level [mJy]")
    pl.savefig(f'{diag_dir}/noise_vs_beam_geomavg.png', bbox_inches='tight')


def estimate_global_rms(data, nbeams=3, threshold=2.0):
    """
    Estimate global continuum RMS using iterative masking.
    """
    fh = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits')
    imgdata = fh[0].data[1557:2557, 8500:9500]
    ww = WCS(fh[0].header)
    pixscale = ww.proj_plane_pixel_area()**0.5
    beam = radio_beam.Beam.from_fits_header(fh[0].header)

    kernelwidth = ((beam.major * nbeams) / pixscale).decompose()

    rms1 = rms_map(imgdata, kernel=Gaussian2DKernel(kernelwidth))
    datacopy = copy.copy(imgdata)
    rms_ = copy.copy(rms1)

    for _ in range(50):  # max iterations
        detections = (datacopy / rms_) > threshold
        if detections.sum() == 0:
            break
        datacopy[detections] = np.nan
        rms_ = rms_map(datacopy, kernel=Gaussian2DKernel(kernelwidth))

    pl.figure()
    pl.imshow(rms_)
    pl.title(f"Threshold = {threshold} sigma")
    pl.colorbar()
    pl.savefig(f'{diag_dir}/global_rms_estimation.png', bbox_inches='tight')


def plot_beamsize_per_field(data):
    high = np.array([('spw33_35' in fn) and ('v1' not in fn) for fn in data['filename']])
    low = np.array([('spw25_27.' in fn) and ('v1' not in fn) for fn in data['filename']])
    agg = np.array([('spw25_27_29_31_33_35' in fn) and ('v1' not in fn) for fn in data['filename']])
    pl.figure(figsize=(12,12))
    toplot = {reg: [row['bmaj']
                    for subset in (agg, low, high)
                    for row in data[(data['region'] == reg) & subset &
                                    (np.char.find(data['filename'], 'preQA3')== -1) &
                                    (np.char.find(data['filename'], 'X15a0_X178')== -1) &
                                    (np.char.find(data['filename'], 'X15a0_Xd6')== -1) &
                                    (np.char.find(data['filename'], 'X15a0_Xe8')== -1) &
                                    (np.char.find(data['filename'], 'obsolete')== -1) &
                                    (np.char.find(data['filename'], 'fits')== -1) &
                                    # either X184 has lower_plus_upper, or it's not X184
                                    (((np.char.find(data['filename'], 'X15a0_X184') >= 0) &
                                      (np.char.find(data['filename'], 'lower_plus_upper') >= 0)) |
                                     ((np.char.find(data['filename'], 'X15a0_X184') == -1))) &
                                    data['pbcor']]]
            for reg in np.unique(data['region'])}
    
    toplot['am'] = [data['bmaj'][(data['suffix'] == 'lower_plus_upper') &
                                 (data['region'] == 'am') &
                                 (data['spws'] == spws)][0]
                    for spws in ('25,27', '25,27,29,31,33,35', '33,35')]
    
    for ii, regname in enumerate(toplot):
        if len(toplot[regname]) > 0: #TODO: this is a hack b/c am is missing
            #print(regname, toplot[regname])
            # sort from small to large to avoid overlap back-and-forth scribbles
            pl.plot(sorted(toplot[regname]), [len(toplot) - ii]*3, color='k', linestyle='--', alpha=0.5, linewidth=0.5)
            for xv, color, label in zip(sorted(toplot[regname]), ['blue', 'red', 'orange',], ('33,35', 'aggregate', '25,27',)):
                pl.plot(xv, len(toplot) - ii, 's', color=color, label=label)
        else:
            print(f"Warning: {regname} has no data")
            raise ValueError(f"ERROR: {regname} has no data")

    pl.legend(loc='best', handles=pl.gca().lines[-3:])
    # "a" should be at the top (len - ii ensures that)
    pl.yticks(np.arange(len(toplot)) + 1,
              np.array(list(toplot.keys()))[::-1], fontsize=14)
    pl.ylim(0.5, len(toplot)+0.5)
    pl.xlabel("Angular Size [\"]")
    pl.savefig(f'{diag_dir}/beamsize_per_field.png', bbox_inches='tight')
    pl.savefig(f'{diag_dir}/beamsize_per_field.pdf', bbox_inches='tight')


def plot_spatial_distribution_of_beams(data):
    fig = pl.figure()
    ax = pl.gca()
    scale = 32
    notfits = (np.char.find(data['filename'], 'fits') == -1)
    TM = data['array'] == '12M'
    high = np.array([('spw33_35' in fn) and ('v1' not in fn) for fn in data['filename']])
    low = np.array([('spw25_27.' in fn) and ('v1' not in fn) for fn in data['filename']])
    agg = np.array([('spw25_27_29_31_33_35' in fn) and ('v1' not in fn) for fn in data['filename']])

    for ii, field in enumerate(np.unique(data['region'][TM])):
        fld = data['region'] == field
        sel = TM & fld & agg & notfits
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if sel.sum() == 0:
                sel = TM & fld & agg & ~notfits
                try:
                    cube = SpectralCube.read(data['filename'][sel][0])
                except:
                    cube = Slice.from_hdu(fits.open(data['filename'][sel][0]))
            else:
                cube = SpectralCube.read(data['filename'][sel][0], format='casa_image')
            ww = cube.wcs.celestial
        try:
            xy = ww.pixel_to_world(cube.shape[2]//2, cube.shape[1]//2).galactic
        except:
            xy = ww.pixel_to_world(cube.shape[1]//2, cube.shape[0]//2).galactic
        ell = pl.matplotlib.patches.Ellipse(
                            xy=(xy.l.wrap_at(180*u.deg).value, xy.b.value),
                            width=data['bmaj'][sel][0] / scale,
                            height=data['bmin'][sel][0] / scale,
                            angle=data['bpa'][sel][0])
        ell.set_edgecolor('k')
        ell.set_facecolor('none')
        ax.text(xy.l.wrap_at(180*u.deg).value, xy.b.value, field, ha='center', va='center', size=12)
        ax.add_patch(ell)

    ell = pl.matplotlib.patches.Ellipse((0.3, -0.18), width=1/scale, height=1/scale, facecolor='none', edgecolor='b')
    ax.add_patch(ell)
    ell = pl.matplotlib.patches.Ellipse((0.3, -0.18), width=3/scale, height=3/scale, facecolor='none', edgecolor='b')
    ax.add_patch(ell)
    ax.text(0.25, -0.25, '1", 3"', ha='left', va='center', color='b')

    ax.axis([0.95, -0.61, -0.3, 0.22])
    ax.set_aspect(1)
    ax.set_xlabel(r"Galactic Longitude ($^\circ$)")
    ax.set_ylabel(r"Galactic Latitude ($^\circ$)")

    fig.savefig(f'{diag_dir}/spatial_distribution_of_beams.pdf',
                bbox_inches='tight') 


def plot_rms_histogram(data):
    pl.figure(figsize=(10, 8))
    rmsmap = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_maskedrms_mosaic.fits')
    beam = radio_beam.Beam.from_fits_header(rmsmap[0].header)
    h1, l1, p1 = pl.hist(rmsmap[0].data[np.isfinite(rmsmap[0].data)]*1e3, bins=np.geomspace(2e-5, 2e-3)*1e3, histtype='step', label=f'Full {beam.major.to(u.arcsec):0.1f}');
    rmsmap = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_maskedrms_mosaic.fits')
    beam = radio_beam.Beam.from_fits_header(rmsmap[0].header)
    h2, l2, p2 = pl.hist(rmsmap[0].data[np.isfinite(rmsmap[0].data)]*1e3, bins=np.geomspace(2e-5, 2e-3)*1e3, histtype='step', label=f'SPW 33+35 {beam.major.to(u.arcsec):0.1f}');
    rmsmap2527 = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_maskedrms_mosaic.fits')
    beam = radio_beam.Beam.from_fits_header(rmsmap2527[0].header)
    h3, l3, p3 = pl.hist(rmsmap2527[0].data[np.isfinite(rmsmap2527[0].data)]*1e3, bins=np.geomspace(2e-5, 2e-3)*1e3, histtype='step', label=f'SPW 25+27 {beam.major.to(u.arcsec):0.1f}');
    pl.semilogx();
    pl.xlabel("RMS Flux Density [mJy beam$^{-1}$]")
    pl.ylabel("Number of 0.2x0.2\" Pixels")
    pl.legend(loc='center right', fancybox=True, framealpha=0.5)
    pl.savefig(f'{diag_dir}/RMS_Histogram_comparison_circular.png', bbox_inches='tight')

    return p1, p2, p3


def plot_rms_histogram_and_cdf(data):
    p1, p2, p3 = plot_rms_histogram(data)
    # CDF
    ax2 = pl.twinx()
    rmsmap = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_maskedrms_mosaic.fits')
    rmssrt = np.sort(rmsmap[0].data[np.isfinite(rmsmap[0].data)])*1e3
    print(f'rmssrt[0] = {rmssrt[0]}.  rmssrt[-1] = {rmssrt[-1]}')
    ax2.plot(rmssrt, np.arange(rmssrt.size)/rmssrt.size, color=p1[0].get_edgecolor(), linestyle='--')

    rmsmap = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_maskedrms_mosaic.fits')
    rmssrt = np.sort(rmsmap[0].data[np.isfinite(rmsmap[0].data)])*1e3
    ax2.plot(rmssrt, np.arange(rmssrt.size)/rmssrt.size,  color=p2[0].get_edgecolor(), linestyle=':')

    rmsmap = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_maskedrms_mosaic.fits')
    rmssrt = np.sort(rmsmap[0].data[np.isfinite(rmsmap[0].data)])*1e3
    ax2.plot(rmssrt, np.arange(rmssrt.size)/rmssrt.size,  color=p3[0].get_edgecolor(), linestyle=':')

    ax2.set_xlim(2e-5*1e3, 2e-3*1e3)
    ax2.set_ylim(0,1)
    ax2.set_ylabel("Fraction of pixels (CDF)")
    #pl.legend(loc='center right')
    pl.savefig(f'{diag_dir}/RMS_Histogram_comparison_circular_withCDF.png', bbox_inches='tight')
    pl.savefig(f'{diag_dir}/RMS_Histogram_comparison_circular_withCDF.pdf', bbox_inches='tight')


def main():
    data = load_and_filter_data()
    plot_beam_size_histograms(data)
    plot_noise_vs_beam(data)
    estimate_global_rms(data)
    plot_beamsize_per_field(data)
    plot_spatial_distribution_of_beams(data)
    plot_rms_histogram_and_cdf(data)

    return data


if __name__ == "__main__":
    data = main()