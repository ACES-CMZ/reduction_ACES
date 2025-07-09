import os
import itertools
import glob
import numpy as np
import matplotlib.pyplot as pl
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

from aces import conf
basepath = conf.basepath

def main():
    contnames = {'continuum': '/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits',
                'continuum_feathered': '/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'}

    for mol1, mol2 in itertools.combinations(['continuum', 'continuum_feathered', "CH3CHO", "CS21", "H13CN", "H13COp", "H40a", "HC15N", "HC3N", "HCOP", "HCOP_mopra", "HCOP_noTP", "HN13C", "HNCO_7m12mTP", "NSplus", "SiO21", "SO21", "SO32",], 2):
        print(f'{mol1} vs {mol2}')
        fn1 = f'{basepath}/mosaics/cubes/moments/{mol1}_CubeMosaic_masked_hlsig_dilated_mom0.fits'
        fn2 = f'{basepath}/mosaics/cubes/moments/{mol2}_CubeMosaic_masked_hlsig_dilated_mom0.fits'

        if mol1 in contnames:
            fn1 = contnames[mol1]
        if mol2 in contnames:
            fn2 = contnames[mol2]

        data1 = fits.getdata(fn1)
        data2 = fits.getdata(fn2)

        fig, ax = pl.subplots(1, 1)
        ax.scatter(data1, data2, s=1, alpha=0.5)
        ax.set_xlabel(f'{mol1} [K km s$^{{-1}}$]')
        ax.set_ylabel(f'{mol2} [K km s$^{{-1}}$]')
        if mol1 in contnames:
            ax.set_xlabel(f'{mol1} [K]')
        if mol2 in contnames:
            ax.set_ylabel(f'{mol2} [K]')

        fig.savefig(f'{basepath}/diagnostic_plots/correlations/{mol1}_{mol2}_correlation.png', bbox_inches='tight', dpi=150)

        pl.close('all')

        fig, ax = pl.subplots(1, 1)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.scatter(data1, data2, s=0.5, alpha=0.05)
        ax.set_xlabel(f'{mol1} [K km s$^{{-1}}$]')
        ax.set_ylabel(f'{mol2} [K km s$^{{-1}}$]')
        if mol1 in contnames:
            ax.set_xlabel(f'{mol1} [K]')
        if mol2 in contnames:
            ax.set_ylabel(f'{mol2} [K]')

        fig.savefig(f'{basepath}/diagnostic_plots/correlations/{mol1}_{mol2}_correlation_log.png', bbox_inches='tight', dpi=150)

        print(f'{mol1} vs {mol2} done')
        pl.close('all')


if __name__ == '__main__':
    main()