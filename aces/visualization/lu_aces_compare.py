#!/usr/bin/env python
"""
Script to compare ALMA data from Xing Lu's projects with ACES data
"""
import os
import sys
import numpy as np
import pylab as pl
from astropy import units as u
from aces.visualization import compare_to_other

from aces.visualization import figure_configuration
from aces.visualization.figure_configuration import mymap, mymap2, distance

from aces import conf


def main():
    """
    Main function to execute the comparisons
    """
    pl.rcParams['font.size'] = 14

    # Calculate effective frequency as in the notebook
    nueff = (86.71 * 0.9375 + 88.29 * 0.9375 + 98.10 * 1.875 + 100 * 1.875) / (1.875 * 2 + 0.9375 * 2)
    print(f"Using effective frequency: {nueff} GHz")

    # Common parameters
    common_kwargs = {
        'wspace': 0.01,
        'hspace': 0.03,
        'norm_kwargs': {'vmax': 1.0, 'vmin': -0.02, 'stretch': 'asinh'},
    }

    # Compare 20kms Cloud
    compare_to_other.compare_to_other_data(
        '/orange/adamginsburg/ACES/ancillary/data/XingLu/C20kms_band3_cont_rebin30_TN2only.image.tt0.pbcor.fits',
        '20kmsCloud',
        nueff*u.GHz,
        figsize=(9, 10),
        **common_kwargs
    )

    # Save the figure
    pl.savefig('20kmsCloud_comparison.pdf', bbox_inches='tight')

    # Compare Cloud E
    compare_to_other.compare_to_other_data(
        '/orange/adamginsburg/ACES/ancillary/data/XingLu/SgrB1off_band3_cont_rebin30_TM2only.image.tt0.pbcor.fits',
        'Cloud E',
        nueff*u.GHz,
        figsize=(12, 10),
        **common_kwargs
    )

    # Save the figure
    pl.savefig('CloudE_comparison.pdf', bbox_inches='tight')

    # Compare SgrC
    compare_to_other.compare_to_other_data(
        '/orange/adamginsburg/ACES/ancillary/data/XingLu/SgrC_band3_cont_rebin30_TM2only.image.tt0.pbcor.fits',
        'SgrC',
        nueff*u.GHz,
        figsize=(12, 10),
        **common_kwargs
    )

    # Save the figure
    pl.savefig('SgrC_comparison.pdf', bbox_inches='tight')


if __name__ == "__main__":
    main()
