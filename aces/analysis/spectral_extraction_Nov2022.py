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

    product_dict = {'SgrC_coreMM1': f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X172/",
                    'SgrC_coreMM2': f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X172/",
                    "BrickMaserCore": f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X190/",
                    "QuakeClump1": f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X1a8/",
                    "QuakeClump2": f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X1a8/",
                    }

    regions_dict = {'SgrC_coreMM1': regions.CircleSkyRegion(coordinates.SkyCoord('17:44:40.60', '-29:28:16.0', frame='icrs', unit=(u.h, u.deg)), radius=0.75 * u.arcsec),
                    'SgrC_coreMM2': regions.CircleSkyRegion(coordinates.SkyCoord('17:44:40.21', '-29:28:12.5', frame='icrs', unit=(u.h, u.deg)), radius=0.75 * u.arcsec),
                    "BrickMaserCore": regions.CircleSkyRegion(coordinates.SkyCoord(0.26087, 0.01616, frame='galactic', unit=(u.deg, u.deg)), radius=0.75 * u.arcsec),
                    "QuakeClump1": regions.CircleSkyRegion(coordinates.SkyCoord('17:46:19.2287619119', '-28:38:29.2978617171', unit=(u.hourangle, u.deg)), radius=0.75 * u.arcsec),
                    "QuakeClump2": regions.CircleSkyRegion(coordinates.SkyCoord('17:46:19.5945777603', '-28:38:09.9402161687', unit=(u.hourangle, u.deg)), radius=0.75 * u.arcsec),
                    }

    for regname in product_dict:
        print()
        print(regname, flush=True)
        product_dir = product_dict[regname]
        reg = regions_dict[regname]

        cubefns = glob.glob(f'{product_dir}/calibrated/working/*Sgr_A_star*.cube.I.iter1*.image.pbcor.fits')

        for fn in cubefns:
            cube = SpectralCube.read(fn)
            print(cube, flush=True)
            cube = cube.subcube_from_regions([reg])
            print(cube, flush=True)
            avg = cube.mean(axis=(1, 2))
            hdu = avg.hdu
            y, x = cube.wcs.celestial.world_to_pixel(reg.center)

            slc = getslice()[int(floor(x)):int(ceil(x)), int(floor(y)):int(ceil(y)), :]
            ww = wcs_utils.slice_wcs(cube.wcs, slc, numpy_order=False)
            hdu.header.update(ww.to_header())
            hdu.data = hdu.data[:, None, None]

            outfn = f"{basepath}/wp3_spectra_Feb2022/{regname}_average_" + fn.split("/")[-1]
            print(outfn, flush=True)
            hdu.writeto(outfn, overwrite=True)
