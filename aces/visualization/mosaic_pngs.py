import time
from spectral_cube import SpectralCube
import os

from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u

from aces import conf
from aces.imaging.make_mosaic import makepng

basepath = conf.basepath

cubepath = f'{basepath}/mosaics/cubes/'
mompath = f'{basepath}/mosaics/cubes/moments/'


def pngify(mol, suf='mom0', **kwargs):
    fh = fits.open(f'{mompath}/{mol}_CubeMosaic_{suf}.fits')
    data = fh[0].data
    wcs = WCS(fh[0].header)

    # DEBUG
    print(f"{mol} {suf} Pixel of 0,0: {wcs.world_to_pixel(SkyCoord(0*u.deg, 0*u.deg, frame='galactic'))}, shape: {data.shape}, crpix={wcs.wcs.crpix}, crval={wcs.wcs.crval}")

    makepng(data=data, wcs=wcs, imfn=f"{mompath}/{mol}_CubeMosaic_{suf}.png",
            **kwargs)


if __name__ == "__main__":
    pngify('SiO21', 'mom0', stretch='sqrt', min_cut=-0.1, max_percent=99.95)
    pngify('SiO21', 'max', stretch='sqrt', min_cut=0, max_percent=99.5)
    pngify('CS21', 'mom0', stretch='sqrt', min_cut=-0.1, max_percent=99.5)
    pngify('CS21', 'max', stretch='sqrt', min_cut=0, max_percent=99.5)
    pngify('HNCO_7m12mTP', 'mom0', stretch='sqrt', min_cut=-0.1, max_percent=99.5)
    pngify('HNCO_7m12mTP', 'max', stretch='sqrt', min_cut=0, max_percent=99.5)
    pngify('HNCO', 'mom0', stretch='sqrt', min_cut=-0.1, max_percent=99.5)
    pngify('HNCO', 'max', stretch='sqrt', min_cut=0, max_percent=99.5)
    pngify('HCOP', 'mom0', stretch='sqrt', min_cut=-0.1, max_percent=99.5)
    pngify('HCOP', 'max', stretch='sqrt', min_cut=0, max_percent=99.5)
    pngify('HC3N', 'mom0', stretch='sqrt', min_cut=-0.1, max_percent=99.5)
    pngify('HC3N', 'max', stretch='sqrt', min_cut=0, max_percent=99.5)