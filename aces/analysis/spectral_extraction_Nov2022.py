import regions
from numpy import floor, ceil
from astropy import coordinates
from astropy import units as u
import glob
from spectral_cube import SpectralCube
from spectral_cube import wcs_utils

from .. import conf
basepath = conf.basepath


class getslice(object):
    def __getitem__(self, x):
        return x


if __name__ == "__main__":

    product_dict = {'SgrC_coreMM1': f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001/X15a0/X172/",
                    'SgrC_coreMM2': f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001/X15a0/X172/"
                    }

    regions_dict = {'SgrC_coreMM1': regions.CircleSkyRegion(coordinates.SkyCoord('17:44:40.60', '-29:28:16.0', frame='icrs', unit=(u.h, u.deg)), radius=0.75 * u.arcsec),
                    'SgrC_coreMM2': regions.CircleSkyRegion(coordinates.SkyCoord('17:44:40.21', '-29:28:12.5', frame='icrs', unit=(u.h, u.deg)), radius=0.75 * u.arcsec),
                    }

    for regname in product_dict:
        print(regname)
        product_dir = product_dict[regname]
        reg = regions_dict[regname]

        cubefns = glob.glob(f'{product_dir}/product/*Sgr_A_star*.cube.I.pbcor.fits')

        for fn in cubefns:
            cube = SpectralCube.read(fn).subcube_from_regions([reg])
            print(cube)
            avg = cube.mean(axis=(1, 2))
            hdu = avg.hdu
            y, x = cube.wcs.celestial.world_to_pixel(reg.center)

            slc = getslice()[int(floor(x)):int(ceil(x)), int(floor(y)):int(ceil(y)), :]
            ww = wcs_utils.slice_wcs(cube.wcs, slc, numpy_order=False)
            hdu.header.update(ww.to_header())
            hdu.data = hdu.data[:, None, None]

            hdu.writeto(f"{regname}_average_" + fn.split("/")[-1], overwrite=True)
