import numpy as np
import regions
import glob
from spectral_cube.spectral_cube import _regionlist_to_single_region
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import SpectralCube
from spectral_cube.cube_utils import mosaic_cubes
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from aces.imaging.make_mosaic import all_lines as all_lines_

from aces import conf

basepath = conf.basepath

# target_wcs = wcs.WCS(naxis=2)
# target_wcs.wcs.ctype = ['GLON-CAR', 'GLAT-CAR']
# target_wcs.wcs.crval = [0, 0]
# target_wcs.wcs.cunit = ['deg', 'deg']
# target_wcs.wcs.cdelt = [-6.38888888888888e-4, 6.388888888888888e-4]
# target_wcs.wcs.crpix = [2000, 2000]
#
# header = target_wcs.to_header()
# header['NAXIS1'] = 4000
# header['NAXIS2'] = 4000

def all_lines(*args, folder='TP_flattened', **kwargs):
    return all_lines_(*args, folder=folder, **kwargs)


def cube_mosaicing():

    for spw in (17, 19, 21, 23, 25, 27):
        filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw{spw}.cube.I.sd.fits')

        # cube0 = SpectralCube.read(filelist[0])
        # output_wcs = WCS(naxis=3)
        # output_wcs.wcs.ctype = ['GLON-CAR', 'GLAT-CAR', 'FREQ']
        # output_wcs.wcs.crval = [0, 0, cube0.wcs.wcs.crval[2]]
        # output_wcs.wcs.cdelt = cube0.wcs.wcs.cdelt
        # output_wcs.wcs.cunit = cube0.wcs.wcs.cunit
        # output_shape = [int(np.abs(1.5 / cube0.wcs.wcs.cdelt[0])),
        #                 int(np.abs(0.55 / cube0.wcs.wcs.cdelt[1])),
        #                 cube0.shape[0]]
        # # crpix * cdelt + crval = cornerval
        # output_wcs.wcs.crpix = [-0.92 / output_wcs.wcs.cdelt[0],
        #                         -0.30 / output_wcs.wcs.cdelt[1],
        #                         cube0.wcs.wcs.crpix[2]]

        result = mosaic_cubes([SpectralCube.read(fn) for fn in filelist],
                              combine_header_kwargs=dict(frame='galactic',
                                                         spectral_dx_threshold=0.01))
        result.write(f'{basepath}/mosaics/cubes/ACES_TP_spw{spw}_mosaic.fits', overwrite=True)


def noisemaps():
    for fn in glob.glob(f"{basepath}/mosaics/cubes/ACES_TP_spw*mosaic.fits"):
        outname = f"{basepath}/mosaics/cubes/moments/{fn.split('/')[-1].replace('mosaic.fits', 'mosaic_madstd.fits')}"
        cube = SpectralCube.read(fn, use_dask=True)
        mstd = cube.mad_std(axis=0)
        try:
            mstd.write(outname, overwrite=True)
        except Exception as ex:
            print(ex)


def main():

    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw17.cube.I.sd.fits')

    def read_as_2d(fn):
        cube = SpectralCube.read(fn, use_dask=True)
        mx = cube.max(axis=0)
        return mx.hdu

    hdus = [read_as_2d(fn) for fn in filelist]

    target_wcs, shape_out = find_optimal_celestial_wcs(hdus, frame='galactic')

    array, footprint = reproject_and_coadd(hdus,
                                           target_wcs, shape_out=shape_out,
                                           reproject_function=reproject_interp)
    header = target_wcs.to_header()
    fits.PrimaryHDU(data=array, header=header).writeto(f'{basepath}/mosaics/TP_spw17mx_mosaic.fits', overwrite=True)
    # re-load the header so it includes NAXIS
    header = fits.getheader(f'{basepath}/mosaics/TP_spw17mx_mosaic.fits')

    import pylab as pl
    from astropy import visualization
    pl.rc('axes', axisbelow=True)
    pl.matplotlib.use('agg')

    front = 10
    back = -10
    fronter = 20

    fig = pl.figure(figsize=(20, 7))
    ax = fig.add_subplot(111, projection=target_wcs)
    im = ax.imshow(array, norm=visualization.simple_norm(array, stretch='asinh'), zorder=front, cmap='gray')
    pl.colorbar(mappable=im)
    ax.coords[0].set_axislabel('Galactic Longitude')
    ax.coords[1].set_axislabel('Galactic Latitude')
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.coords[0].set_ticks(spacing=0.1 * u.deg)
    ax.coords[0].set_ticklabel(rotation=45, pad=20)

    fig.savefig(f'{basepath}/mosaics/TP_spw17mx_mosaic.png', bbox_inches='tight')

    ax.coords.grid(True, color='black', ls='--', zorder=back)

    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='black', ls=':', zorder=back)
    overlay[0].set_axislabel('Right Ascension (ICRS)')
    overlay[1].set_axislabel('Declination (ICRS)')
    overlay[0].set_major_formatter('hh:mm')
    ax.set_axisbelow(True)
    ax.set_zorder(back)
    fig.savefig(f'{basepath}/mosaics/TP_spw17mx_mosaic_withgrid.png', bbox_inches='tight')

    tbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/SB_naming.md', format='ascii.csv', delimiter='|', data_start=2)

    # create a list of composite regions for labeling purposes
    composites = []
    flagmap = np.zeros(array.shape, dtype='int')
    # loop over the SB_naming table
    for row in tbl:
        indx = row['Proposal ID'][3:]

        # load up regions
        regs = regions.Regions.read(f'{basepath}/reduction_ACES/aces/data/regions/final_cmz{indx}.reg')
        pregs = [reg.to_pixel(target_wcs) for reg in regs]

        composite = _regionlist_to_single_region(pregs)
        composite.meta['label'] = indx
        composites.append(composite)

        for comp in composites:
            cmsk = comp.to_mask()

            slcs_big, slcs_small = cmsk.get_overlap_slices(array.shape)
            try:
                flagmap[slcs_big] += (cmsk.data[slcs_small] * int(comp.meta['label'])) * (flagmap[slcs_big] == 0)
            except ValueError:
                # expected to occur if no overlap
                continue

    ax.contour(flagmap, cmap='prism', levels=np.arange(flagmap.max()) + 0.5, zorder=fronter)

    for ii in np.unique(flagmap):
        if ii > 0:
            fsum = (flagmap == ii).sum()
            cy, cx = ((np.arange(flagmap.shape[0])[:, None] * (flagmap == ii)).sum() / fsum,
                      (np.arange(flagmap.shape[1])[None, :] * (flagmap == ii)).sum() / fsum)
            pl.text(cx, cy, f"{ii}\n{tbl[ii - 1]['Obs ID']}",
                    horizontalalignment='left', verticalalignment='center',
                    color=(1, 0.8, 0.5), transform=ax.get_transform('pixel'),
                    zorder=fronter)

    fig.savefig(f'{basepath}/mosaics/TP_spw17mx_mosaic_withgridandlabels.png', bbox_inches='tight')

    cube_mosaicing()
    noisemaps()

    all_lines(header, array='TP', glob_suffix='cube.I.sd.fits', globdir='product', use_weights=False,
              parallel=False)