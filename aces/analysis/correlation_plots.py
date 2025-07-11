import os
import time
import itertools
import glob
import numpy as np
import matplotlib.pyplot as pl
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import mpl_plot_templates
import regions

from aces import conf
basepath = conf.basepath

def main():
    os.makedirs(f'{basepath}/diagnostic_plots/correlations', exist_ok=True)

    contnames = {'continuum': '/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits',
                'continuum_feathered': '/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'}

    t0 = time.time()

    ww = WCS(fits.getheader(contnames['continuum_feathered']))

    pixreg = regions.Regions.read(f'{basepath}/regions/sgramask.reg')[0].to_pixel(ww).union(
        regions.Regions.read(f'{basepath}/regions/sgrb2mask.reg')[0].to_pixel(ww))
    reg_mask = pixreg.to_mask()
    sgra_or_sgrb2 = (reg_mask.to_image(fits.getdata(contnames['continuum_feathered']).shape) > 0)
    print(sgra_or_sgrb2.shape, flush=True)


    #for mol1, mol2 in itertools.combinations([ "H40a",  "SO21", "SO32",], 2):
    for mol1, mol2 in itertools.combinations(['continuum', 'continuum_feathered', "CH3CHO", "CS21", "H13CN", "H13COp", "H40a", "HC15N", "HC3N", "HCOP", "HCOP_mopra", "HCOP_noTP", "HN13C", "HNCO_7m12mTP", "NSplus", "SiO21", "SO21", "SO32",], 2):
        print(f'{mol1} vs {mol2}.   time={time.time() - t0:0.2f}', flush=True)
        fn1 = f'{basepath}/mosaics/cubes/moments/{mol1}_CubeMosaic_masked_hlsig_dilated_mom0.fits'
        fn2 = f'{basepath}/mosaics/cubes/moments/{mol2}_CubeMosaic_masked_hlsig_dilated_mom0.fits'

        if mol1 in contnames:
            fn1 = contnames[mol1]
        if mol2 in contnames:
            fn2 = contnames[mol2]

        data1 = fits.getdata(fn1)
        data2 = fits.getdata(fn2)
        print(data1.shape, data2.shape)

        to_xcorr = np.isfinite(data1) & np.isfinite(data2)
        corr = np.corrcoef(data1[to_xcorr], data2[to_xcorr])[0, 1]
        assert np.isfinite(corr)

        xmin, xmax = np.nanmin(data1), np.nanmax(data1)
        x5, x95 = np.nanpercentile(data1, [5, 95])
        y5, y95 = np.nanpercentile(data2, [5, 95])

        fig, ax = pl.subplots(1, 1)
        ax.plot([xmin, xmax], [xmin, xmax], color='b', linestyle='--')
        ax.scatter(data1, data2, s=1, alpha=0.5, color='k', label=f'r={corr:.2f}')
        ax.scatter(data1[sgra_or_sgrb2], data2[sgra_or_sgrb2], s=1, alpha=0.5, color='r', label=f'Sgr B2 & Sgr A*')
        ax.set_xlabel(f'{mol1} [K km s$^{{-1}}$]')
        ax.set_ylabel(f'{mol2} [K km s$^{{-1}}$]')
        if mol1 in contnames:
            ax.set_xlabel(f'{mol1} [K]')
        if mol2 in contnames:
            ax.set_ylabel(f'{mol2} [K]')
        ax.legend()

        fig.savefig(f'{basepath}/diagnostic_plots/correlations/{mol1}_{mol2}_correlation.png', bbox_inches='tight', dpi=150)

        pl.close('all')

        bothpos = (data1 > 0) & (data2 > 0) & ~sgra_or_sgrb2
        corr = np.corrcoef(data1[bothpos], data2[bothpos])[0, 1]

        fig, ax = pl.subplots(1, 1)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.hexbin(data1[bothpos], data2[bothpos], xscale='log', yscale='log', mincnt=1)
        ax.plot([xmin, xmax], [xmin, xmax], color='k', linestyle='--', label=f'r={corr:.2f}')
        #ax.scatter(data1, data2, s=0.1, alpha=0.05, color='k')
        # rslt = mpl_plot_templates.adaptive_param_plot(x=data1, y=data2,
        #                                               bins=[np.geomspace(x5, x95, 25), np.geomspace(y5, y95, 25)],
        #                                               threshold=30, marker=',')

        ax.scatter(data1[sgra_or_sgrb2 & bothpos], data2[sgra_or_sgrb2 & bothpos], s=1, alpha=0.5, color='r', label=f'Sgr B2 & Sgr A*')

        ax.set_xlabel(f'{mol1} [K km s$^{{-1}}$]')
        ax.set_ylabel(f'{mol2} [K km s$^{{-1}}$]')
        if mol1 in contnames:
            ax.set_xlabel(f'{mol1} [K]')
        if mol2 in contnames:
            ax.set_ylabel(f'{mol2} [K]')
        ax.legend()

        fig.savefig(f'{basepath}/diagnostic_plots/correlations/{mol1}_{mol2}_correlation_log.png', bbox_inches='tight', dpi=150)

        print(f'{mol1} vs {mol2} done', flush=True)
        pl.close('all')


if __name__ == '__main__':
    main()