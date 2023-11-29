"""
sbatch --qos=astronomy-dept-b --account=astronomy-dept --ntasks=1 --nodes=1 --mem=64gb --time=96:00:00 --job-name=ACES_fixMultiBeam --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/casapy38/bin/python /orange/adamginsburg/ACES/reduction_ACES/aces/hipergator_scripts/fix_multibeam_Aug1_2023.py"
"""
import os
import shutil
import glob
from casatasks import imsmooth, impbcor, exportfits
from casatools import image
ia = image()

images = glob.glob("/orange/adamginsburg/ACES/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/m*/calibrated/working/*.image")

for dotimage in images:
    imagename = os.path.splitext(dotimage)[0]

    os.chdir(os.path.dirname(imagename))

    try:
        psffile = os.path.basename(f'{imagename}.psf')
        ia.open(psffile)
        commonbeam = ia.commonbeam()
        ia.close()
    except Exception as ex:
        print(f"Problem with {psffile}: {ex}")
        ia.open(dotimage)
        commonbeam = ia.commonbeam()
        ia.close()

    ia.open(dotimage)
    rbeam = ia.restoringbeam()
    ia.close()

    # if any beam is not the common beam...
    manybeam = False
    if 'beams' in rbeam:
        for beam in rbeam['beams'].values():
            if beam['*0']['major'] != commonbeam['major'] or beam['*0']['minor'] != commonbeam['minor']:
                print(f"beam={beam['*0']}, commonbeam={commonbeam}")
                manybeam = True
                break

    if manybeam:
        print(f"Found a multi-beam image {imagename} = {dotimage}")
        if os.path.exists(f'{imagename}.image.multibeam'):
            raise ValueError(f"The multi-beam image {dotimage} should have already been corrected")
        shutil.move(f'{imagename}.image', f'{imagename}.image.multibeam')
        shutil.move(f'{imagename}.image.pbcor', f'{imagename}.image.pbcor.multibeam')
        imsmooth(imagename=f'{imagename}.model',
                 outfile=f'{imagename}.convmodel',
                 beam=commonbeam)
        ia.imagecalc(outfile=f'{os.path.basename(imagename)}.image',
                     pixels=f'{os.path.basename(imagename)}.convmodel + {os.path.basename(imagename)}.residual',
                     imagemd=f'{os.path.basename(imagename)}.convmodel',
                     overwrite=True)
        ia.close()
        impbcor(imagename=f'{os.path.basename(imagename)}.image',
                pbimage=f'{os.path.basename(imagename)}.pb',
                outfile=f'{os.path.basename(imagename)}.image.pbcor',)

        exportfits(imagename=f'{os.path.basename(imagename)}.image.pbcor',
                   fitsimage=f'{os.path.basename(imagename)}.image.pbcor.fits',
                   overwrite=True)
        exportfits(imagename=f'{os.path.basename(imagename)}.image',
                   fitsimage=f'{os.path.basename(imagename)}.image.fits',
                   overwrite=True)
