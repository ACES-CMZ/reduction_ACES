import regions
from numpy import floor, ceil
from astropy import coordinates
from astropy import units as u
import glob
from spectral_cube import SpectralCube
from spectral_cube import wcs_utils

from aces import conf
basepath = conf.basepath


class getslice(object):
    def __getitem__(self, x):
        return x


if __name__ == "__main__":

    product_dict = {'brick_h2o_core': f'{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X192/',
                    'SgrB2_G0.69': f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_Xea/",
                    'cloud_e_hotcore': f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X1a4/",
                    'SgrC_hotcore': f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001/X15a0/X174/",
                    'SgrC_MM1_HC_12m': f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001/X15a0/X172/",
                    'brick_maser_core_12m': f'{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X190/',
                    }

    regions_dict = {'brick_h2o_core': regions.CircleSkyRegion(coordinates.SkyCoord('17:46:10.628', '-28:42:17.75', frame='icrs', unit=(u.h, u.deg)), radius=15 * u.arcsec),
                    'SgrB2_G0.69': regions.CircleSkyRegion(coordinates.SkyCoord('17h47m22s', '-28:21:27', frame='fk5', unit=(u.h, u.deg)), radius=15 * u.arcsec),
                    'cloud_e_hotcore': regions.CircleSkyRegion(coordinates.SkyCoord(0.4751733 * u.deg, -0.0096808 * u.deg, frame='galactic'), radius=15 * u.arcsec),
                    'SgrC_hotcore': regions.CircleSkyRegion(coordinates.SkyCoord('17:44:40.58', '-29:28:16.28', frame='icrs', unit=(u.h, u.deg)), radius=15 * u.arcsec),
                    'SgrC_MM1_HC_12m': regions.CircleSkyRegion(coordinates.SkyCoord('17:44:40.58', '-29:28:15.82', frame='icrs', unit=(u.h, u.deg)), radius=3 * u.arcsec),
                    'brick_maser_core_12m': regions.CircleSkyRegion(coordinates.SkyCoord('17:46:10.63', '-28:42:17.77', frame='icrs', unit=(u.h, u.deg)), radius=3 * u.arcsec),
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
