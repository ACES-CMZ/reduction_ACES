import regions
import os
import warnings
from numpy import floor, ceil
from astropy import coordinates
from astropy import units as u
import glob
from spectral_cube import SpectralCube
from spectral_cube import wcs_utils
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import reproject
import scipy.ndimage

from aces import conf
basepath = conf.basepath


class getslice(object):
    def __getitem__(self, x):
        return x


def extract_from_mask(cube, maskhdu, maskid):
    """
    """
    mask_match = maskhdu.data == maskid
    obj_slc = scipy.ndimage.find_objects(mask_match)[0]
    mask_co = mask_match[obj_slc]
    mask_ww = WCS(maskhdu.header)[obj_slc]

    mask_rep, _ = reproject.reproject_interp((mask_co, mask_ww),
                                             cube.wcs.celestial,
                                             shape_out=cube.shape[1:])

    slcs = cube.subcube_slices_from_mask(mask_rep > 0, spatial_only=True)
    scube = cube[slcs]

    spec = scube.mean(axis=(1, 2))

    return spec


if __name__ == "__main__":

    maskfile = fits.open(f"{basepath}/upload/Cont_catalog_stuff/ACES_leaf_mask_3_1_mp179.fits")
    catalog = Table.read(f"{basepath}/upload/Cont_catalog_stuff/aces_catalog_3_1_mp179.fits")

    product_dir = f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/"
    spectrum_dir = f"{basepath}/spectra"

    field_map = fits.open(f'{basepath}/mosaics/12m_continuum_field_number_map.fits')
    fieldmapwcs = WCS(field_map[0].header)
    fieldmapdata = field_map[0].data

    uidtbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/aces_SB_uids.csv')

    for row in catalog:
        center = SkyCoord(row['GLON'], row['GLAT'], frame='galactic', unit=(u.deg, u.deg))
        reg = regions.EllipseSkyRegion(center,
                                       width=row['major_sigma'] * u.arcsec,
                                       height=row['minor_sigma'] * u.arcsec,
                                       angle=row['position_angle'] * u.deg)

        pixcrd = list(map(lambda x: int(x), fieldmapwcs.world_to_pixel(center)))
        field_id = fieldmapdata[pixcrd[1], pixcrd[0]]
        if field_id == 0:
            #print(f"Skipped row {row}: out of bounds")
            continue
        mousid = uidtbl[field_id - 1]['12m MOUS ID']

        # Can make this step more efficient by using the table to map from position to cube ID
        # (right now this does a lot of extra looping & reading and is therefore pretty slow)
        # basically use mosaics/12m_continuum_field_number_map.fits to load reduction_ACES/aces/data/tables/aces_SB_uids.csv
        # to then load the appropriate folder
        cubefns = glob.glob(f'{product_dir}/member.uid___A001_{mousid}/calibrated/working/*Sgr_A_star*.cube.I.iter1.image.pbcor.fits')
        # print(row['index'], field_id, uidtbl[field_id]['12m MOUS ID'])
        # print(cubefns)
        # break
        for cubefn in cubefns:

            # there are two fns to check, but we're only checking the second
            outfn = f"{spectrum_dir}/mp179_source{row['index']}_dendromaskaverage_" + cubefn.split("/")[-1]

            if not os.path.exists(outfn):

                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    cube = SpectralCube.read(cubefn)
                    cube.allow_huge_operations = True  # suppress future warnings
                    ww = cube.wcs.celestial
                    ww._naxis = cube.shape[1:]

                if ww.footprint_contains(center):
                    scube = cube.subcube_from_regions([reg])
                    print(row['index'], cubefn)
                    avg = scube.mean(axis=(1, 2))
                    hdu = avg.hdu

                    y, x = scube.wcs.celestial.world_to_pixel(reg.center)
                    slc = getslice()[int(floor(x)):int(ceil(x)), int(floor(y)):int(ceil(y)), :]
                    ww = wcs_utils.slice_wcs(scube.wcs, slc, numpy_order=False)
                    hdu.header.update(ww.to_header())
                    hdu.data = hdu.data[:, None, None]
                    hdu.header['CATINDX'] = row['index']
                    hdu.header['CATGLON'] = row['GLON']
                    hdu.header['CATGLAT'] = row['GLAT']
                    hdu.header['CATMAJS'] = row['major_sigma']
                    hdu.header['CATMINS'] = row['minor_sigma']
                    hdu.header['CATPA'] = row['position_angle']

                    outfn = f"{spectrum_dir}/mp179_source{row['index']}_ellipseaverage_" + cubefn.split("/")[-1]
                    hdu.writeto(outfn, overwrite=True)

                    try:
                        spec = extract_from_mask(cube, maskfile[0], row['index'] + 1)
                    except Exception as ex:
                        print(f"Failed for cube {cubefn} for id {row['index']} with exception {ex}")
                        continue
                    hdu = spec.hdu
                    hdu.header.update(ww.to_header())
                    hdu.data = hdu.data[:, None, None]
                    hdu.header['CATINDX'] = row['index']
                    hdu.header['CATGLON'] = row['GLON']
                    hdu.header['CATGLAT'] = row['GLAT']
                    hdu.header['CATMAJS'] = row['major_sigma']
                    hdu.header['CATMINS'] = row['minor_sigma']
                    hdu.header['CATPA'] = row['position_angle']

                    outfn = f"{spectrum_dir}/mp179_source{row['index']}_dendromaskaverage_" + cubefn.split("/")[-1]
                    hdu.writeto(outfn, overwrite=True)
                else:
                    print('WTF?', row['index'], cubefn)
