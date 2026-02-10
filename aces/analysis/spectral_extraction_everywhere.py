import regions
import os
import warnings
import numpy as np
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


def extract_background_spectrum(cube, center, source_major, source_minor,
                                inner_scale=1.5, outer_scale=3.0):
    """Extract a background annulus spectrum around a source.
    
    Parameters
    ----------
    cube : SpectralCube
        The spectral cube.
    center : SkyCoord
        Center of the source.
    source_major : float
        Source major axis in arcsec.
    source_minor : float
        Source minor axis in arcsec.
    inner_scale : float
        Scale factor for inner radius (relative to max(major, minor)).
    outer_scale : float
        Scale factor for outer radius (relative to max(major, minor)).
        
    Returns
    -------
    bg_spectrum : Spectrum1D-like
        Background spectrum (mean of annulus).
    """
    # Use circular annulus based on larger axis
    radius = max(source_major, source_minor)
    inner_radius = radius * inner_scale * u.arcsec
    outer_radius = radius * outer_scale * u.arcsec
    
    # Create inner and outer circular regions
    inner_reg = regions.CircleSkyRegion(center, radius=inner_radius)
    outer_reg = regions.CircleSkyRegion(center, radius=outer_radius)
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # Extract outer region
        outer_scube = cube.subcube_from_regions([outer_reg])
        outer_mean = outer_scube.mean(axis=(1, 2))
        
        # Extract inner region
        inner_scube = cube.subcube_from_regions([inner_reg])
        inner_mean = inner_scube.mean(axis=(1, 2))
    
    # Calculate annulus spectrum with proper weighting
    outer_mask = outer_scube.mask.include()
    inner_mask = inner_scube.mask.include()
    n_outer = np.sum(outer_mask, axis=(1, 2)).astype(float)
    n_inner = np.sum(inner_mask, axis=(1, 2)).astype(float)
    n_annulus = n_outer - n_inner
    n_annulus[n_annulus <= 0] = 1.0
    
    bg_data = ((outer_mean.value * n_outer - inner_mean.value * n_inner)
               / n_annulus)
    
    # Create a spectrum with the same structure as outer_mean
    bg_spectrum = outer_mean.copy()
    bg_spectrum._data = bg_data
    
    return bg_spectrum


if __name__ == "__main__":

    #maskfile = fits.open(f"{basepath}/upload/Cont_catalog_stuff/ACES_leaf_mask_3_1_mp179.fits")
    #catalog = Table.read(f"{basepath}/upload/Cont_catalog_stuff/aces_catalog_3_1_mp179.fits")
    catalog = Table.read(f'{basepath}/upload/compact_cont_source_catalog/official_cont_catalog_files/aces_compact_catalog_v0.fits')
    catalog_name_prefix = 'ACEScatalog_v0_20260130'

    product_dir = f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/"
    spectrum_dir = f"{basepath}/spectra"

    field_map = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_field_number_map.fits')
    fieldmapwcs = WCS(field_map[0].header)
    fieldmapdata = field_map[0].data

    uidtbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/aces_SB_uids.csv')

    catalog.add_column([' ' * 60]*len(catalog), name='homecube')

    for row in catalog:
        center = SkyCoord(row['GLON_peak'], row['GLAT_peak'], frame='galactic', unit=(u.deg, u.deg))
        reg = regions.EllipseSkyRegion(center,
                                       width=row['fitted_major'] * u.arcsec,
                                       height=row['fitted_minor'] * u.arcsec,
                                       angle=row['pa'] * u.deg)
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

            # there are two fns to check, but we're only checking the first
            outfn = f"{spectrum_dir}/{catalog_name_prefix}_source{row['index']}_ellipseaverage_" + cubefn.split("/")[-1]

            if os.path.exists(outfn):
                try:
                    # test that the file is openable
                    fh = fits.open(outfn)
                    fh.close()
                    print(f"Skipping existing {outfn}", flush=True)
                    row['homecube'] = os.path.basename(cubefn)
                    continue
                except OSError:
                    print(f"Re-extracting existing {outfn} because it failed to open", flush=True)

            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                cube = SpectralCube.read(cubefn)
                cube.allow_huge_operations = True  # suppress future warnings
                ww = cube.wcs.celestial
                ww._naxis = cube.shape[1:]

            if ww.footprint_contains(center):
                print(row['index'], cubefn)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    scube = cube.subcube_from_regions([reg])
                    avg = scube.mean(axis=(1, 2))
                hdu = avg.hdu

                y, x = scube.wcs.celestial.world_to_pixel(reg.center)
                slc = getslice()[int(floor(x)):int(ceil(x)), int(floor(y)):int(ceil(y)), :]
                ww = wcs_utils.slice_wcs(scube.wcs, slc, numpy_order=False)
                hdu.header.update(ww.to_header())
                hdu.data = hdu.data[:, None, None]
                hdu.header['CATINDX'] = row['index']
                hdu.header['CATGLON'] = row['GLON_peak']
                hdu.header['CATGLAT'] = row['GLAT_peak']
                hdu.header['CATMAJS'] = row['fitted_major']
                hdu.header['CATMINS'] = row['fitted_minor']
                hdu.header['CATPA'] = row['pa']

                outfn = f"{spectrum_dir}/{catalog_name_prefix}_source{row['index']}_ellipseaverage_" + cubefn.split("/")[-1]
                hdu.writeto(outfn, overwrite=True)
                print(outfn, flush=True)
                row['homecube'] = cubefn
                
                # Extract and save background spectrum
                bg_outfn = outfn.replace('.fits', '_background.fits')
                if not os.path.exists(bg_outfn):
                    bg_spectrum = extract_background_spectrum(
                        cube, center,
                        row['fitted_major'], row['fitted_minor']
                    )
                    bg_hdu = bg_spectrum.hdu
                    bg_hdu.header.update(ww.to_header())
                    bg_hdu.data = bg_hdu.data[:, None, None]
                    bg_hdu.header['CATINDX'] = row['index']
                    bg_hdu.header['CATGLON'] = row['GLON_peak']
                    bg_hdu.header['CATGLAT'] = row['GLAT_peak']
                    bg_hdu.header['CATMAJS'] = row['fitted_major']
                    bg_hdu.header['CATMINS'] = row['fitted_minor']
                    bg_hdu.header['CATPA'] = row['pa']
                    bg_hdu.header['BKGTYPE'] = 'annulus'
                    bg_hdu.header['BKGINNER'] = (1.5, 'Inner radius scale factor')
                    bg_hdu.header['BKGOUTER'] = (3.0, 'Outer radius scale factor')
                    bg_hdu.writeto(bg_outfn, overwrite=True)
                    print(f"  Background: {bg_outfn}", flush=True)

                # old version with masks ("final" catalog doesn't have masks)
                # try:
                #     spec = extract_from_mask(cube, maskfile[0], row['index'] + 1)
                # except Exception as ex:
                #     print(f"Failed for cube {cubefn} for id {row['index']} with exception {ex}")
                #     continue
                # hdu = spec.hdu
                # hdu.header.update(ww.to_header())
                # hdu.data = hdu.data[:, None, None]
                # hdu.header['CATINDX'] = row['index']
                # hdu.header['CATGLON'] = row['GLON']
                # hdu.header['CATGLAT'] = row['GLAT']
                # hdu.header['CATMAJS'] = row['major_sigma']
                # hdu.header['CATMINS'] = row['minor_sigma']
                # hdu.header['CATPA'] = row['position_angle']

                # outfn = f"{spectrum_dir}/{catalog_name_prefix}_source{row['index']}_dendromaskaverage_" + cubefn.split("/")[-1]
                # hdu.writeto(outfn, overwrite=True)
        if row['homecube']:
            print(f"Source {row['index']} used cube {row['homecube']}", flush=True)
        else:
            print(f"ERROR: Source {row['index']} was not in any cube", flush=True)

    print(f"There were {len(catalog)} total sources in the catalog, of which {np.sum(catalog['homecube'] != ' ' * 60)} were not in a cube", flush=True)
