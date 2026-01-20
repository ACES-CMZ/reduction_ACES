"""
Link files from /orange/adamginsburg/ACES/mosaics/cubes/ to /orange/adamginsburg/ACES/products_for_ALMA/group.uid___A001_X1590_X30a9.lp_slongmore.galactic
"""
import os
import datetime

molnamelist = ['HCOP_mopra', 'HNCO_7m12mTP', 'SiO21', 'SO21', 'H13COp', 'H13CN', 'HN13C', 'HC15N', 'CS21', 'H40a', 'SO32', 'HC3N', 'CH3CHO', 'NSplus']

srcpath = '/orange/adamginsburg/ACES/mosaics/cubes/'
destpath = '/orange/adamginsburg/ACES/products_for_ALMA/group.uid___A001_X1590_X30a9.lp_slongmore.galactic/'

suffixes_to_targets = {'': 'pbcor',
                       '_downsampled9': 'downsampled_spatially.pbcor',
                       '_downsampled9_spectrally': 'downsampled_spectrally.pbcor'}

timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
# Track which files were linked in case we need to roll back
with open(f'/blue/adamginsburg/adamginsburg/ACES/logs/files_that_were_linked_{timestamp}.txt', 'w') as logfile:

    for suffix, target_suffix in suffixes_to_targets.items():
        for molname in molnamelist:
            src = f'{srcpath}/{molname}_CubeMosaic{suffix}.fits'
            dest = f'{destpath}/group.uid___A001_X1590_X30a9.lp_slongmore.cmz_mosaic.12m7mTP.{molname}.cube.{target_suffix}.fits'

            try:
                os.link(src, dest)
                logfile.write(dest + '\n')
                print(f'Linked: {os.path.basename(src)}->{os.path.basename(dest)}')
            except FileExistsError:
                print(f'File exists, skipping: {os.path.basename(dest)}')