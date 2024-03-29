"""
sbatch --qos=astronomy-dept-b --account=astronomy-dept --ntasks=1 --nodes=1 --mem=64gb --time=96:00:00 --job-name=ACES_fixMultiBeam --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/casapy38/bin/python /orange/adamginsburg/ACES/reduction_ACES/aces/hipergator_scripts/fix_multibeam_Aug1_2023.py"
"""
import os
import numpy as np
import shutil
import glob
from casatasks import imsmooth, impbcor, exportfits
from casatools import image
ia = image()

images = glob.glob("/orange/adamginsburg/ACES/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/m*/calibrated/working/*cube*.image")
#images = glob.glob("/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X1a2/calibrated/working/*spw33.cube*psf")


# if we use .image and the code below fails halfway, it won't try to recover b/c it renames .image to .image.multibeam
# so we can use .residual instead, maybe (or .model)
images = glob.glob("/orange/adamginsburg/ACES/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/m*/calibrated/working/*cube*.residual")

for dotsomething in images:
    print(f"Working on {dotsomething} (fixing multibeam if needed)", flush=True)
    imagename = os.path.splitext(dotsomething)[0]
    dotimage = f"{imagename}.image"

    os.chdir(os.path.dirname(imagename))

    try:
        psffile = os.path.basename(f'{imagename}.psf')
        ia.open(psffile)
        ia.commonbeam()  # this checks whether there is *any* beam, and fails to .image otherwise
        used = 'psf'
    except Exception as ex:
        print(f"Problem with {psffile}: {ex}", flush=True)
        ia.open(dotimage)
        used = 'image'

    beams = ia.restoringbeam()
    if 'beams' in beams:
        bmaj = np.array([bm['*0']['major']['value'] for bm in beams['beams'].values()])
        bad_beams = (bmaj > 1.1 * np.median(bmaj)) | (bmaj < 0.9 * np.median(bmaj))
        if any(bad_beams):
            print(f"Found {bad_beams.sum()} bad beams in {imagename}", flush=True)
            avbeam = np.median(bmaj[~bad_beams])
            for ii, (beam, isbad) in enumerate(zip(beams['beams'].values(), bad_beams)):
                if isbad:
                    beam['*0']['major']['value'] = avbeam
                    beam['*0']['minor']['value'] = avbeam
                    beam['*0']['positionangle']['value'] = 0
                ia.setrestoringbeam(beam=beam['*0'], channel=ii)

    commonbeam = ia.commonbeam()
    ia.close()

    try:
        ia.open(dotimage)
        rbeam = ia.restoringbeam()
        ia.close()
        noimage = False
    except Exception as ex:  # noqa
        rbeam = beams
        noimage = True

    # if any beam is not the common beam...
    manybeam = False
    if 'beams' in rbeam:
        for beam in rbeam['beams'].values():
            if beam['*0']['major'] != commonbeam['major'] or beam['*0']['minor'] != commonbeam['minor']:
                print(f"beam={beam['*0']}, commonbeam={commonbeam}", flush=True)
                manybeam = True
                break

    if manybeam or noimage:
        if noimage:
            print(f"Found no image: {dotimage}", flush=True)
        else:
            print(f"Found a multi-beam image {imagename} = {dotimage}", flush=True)
            if os.path.exists(f'{imagename}.image.multibeam'):
                print(f"The multi-beam image {dotimage} should have already been corrected.  Skipping the move", flush=True)
                #raise ValueError(f"The multi-beam image {dotimage} should have already been corrected")
            else:
                shutil.move(f'{imagename}.image', f'{imagename}.image.multibeam')
                shutil.move(f'{imagename}.image.pbcor', f'{imagename}.image.pbcor.multibeam')

        if not os.path.exists(f'{imagename}.model'):
            print(f"MISSING MODEL for {imagename}")
            continue
        imsmooth(imagename=f'{imagename}.model',
                 outfile=f'{imagename}.convmodel',
                 beam=commonbeam,
                 overwrite=True
                 )
        ia.imagecalc(outfile=f'{os.path.basename(imagename)}.image',
                     pixels=f'{os.path.basename(imagename)}.convmodel + {os.path.basename(imagename)}.residual',
                     imagemd=f'{os.path.basename(imagename)}.convmodel',
                     overwrite=True)
        ia.close()
        impbcor(imagename=f'{os.path.basename(imagename)}.image',
                pbimage=f'{os.path.basename(imagename)}.pb',
                outfile=f'{os.path.basename(imagename)}.image.pbcor',
                overwrite=True
                )

        exportfits(imagename=f'{os.path.basename(imagename)}.image.pbcor',
                   fitsimage=f'{os.path.basename(imagename)}.image.pbcor.fits',
                   overwrite=True)
        exportfits(imagename=f'{os.path.basename(imagename)}.image',
                   fitsimage=f'{os.path.basename(imagename)}.image.fits',
                   overwrite=True)

    # separate step: make sure pbcor exists
    # TODO: check the .fits file also
    pbcimg = f'{os.path.basename(imagename)}.image.pbcor'
    try:
        # check whether there's a single beam
        ia.open(pbcimg)
        rbeam_pbc = ia.restoringbeam()
        ia.close()
        nopbcimage = False
    except Exception as ex:  # noqa
        nopbcimage = True

    # if any beam is not the common beam...
    manybeam_pbc = False
    if 'beams' in rbeam_pbc:
        for beam in rbeam_pbc['beams'].values():
            if beam['*0']['major'] != commonbeam['major'] or beam['*0']['minor'] != commonbeam['minor']:
                print(f"beam={beam['*0']}, commonbeam={commonbeam}", flush=True)
                manybeam_pbc = True
                break

    if nopbcimage or manybeam_pbc or not os.path.exists(pbcimg):
        print(f"pbcorrecting {imagename}", flush=True)
        impbcor(imagename=f'{os.path.basename(imagename)}.image',
                pbimage=f'{os.path.basename(imagename)}.pb',
                outfile=f'{os.path.basename(imagename)}.image.pbcor',
                overwrite=True
                )
        exportfits(imagename=f'{os.path.basename(imagename)}.image.pbcor',
                   fitsimage=f'{os.path.basename(imagename)}.image.pbcor.fits',
                   overwrite=True)