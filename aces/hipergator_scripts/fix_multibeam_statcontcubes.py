from spectral_cube import SpectralCube
from spectral_cube.utils import NoBeamError
from astropy.io import fits
import glob
files = glob.glob('/orange/adamginsburg/ACES//rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_*/calibrated/working/*.statcont.contsub.fits')

def has_beam(cube):
    try:
        cube.beam
        return True
    except (AttributeError, NoBeamError):
        return False

for fn in files:
    try:
        cube = SpectralCube.read(fn)
    except Exception as ex:
        print(f"Cube {fn} failed with error {ex}")
    if hasattr(cube, 'beams'):
        print(f"{fn} is a VaryingResolutionCube!")
    elif not has_beam(cube):
        print(fn)
        try:
            bcube = SpectralCube.read(fn.replace(".statcont.contsub",""))
        except Exception as ex:
            print(f'Cube {fn.replace(".statcont.contsub","")} failed with error {ex}')
            continue
        if has_beam(bcube):
            beam = bcube.beam
            with fits.open(fn, mode='update') as fh:
                for kw in bcube.beam.to_header_keywords():
                    fh[0].header[kw] = bcube.beam.to_header_keywords()[kw]
                fh.flush()
                print(f"Updated kws in {fn}.  bmaj={fh[0].header['BMAJ']}")
        elif hasattr(bcube, 'beams'):
            print(f"{fn} is based on a VaryingResolutionCube")
