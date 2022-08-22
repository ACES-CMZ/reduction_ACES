import numpy as np
import regions
import radio_beam
import spectral_cube
from spectral_cube.spectral_cube import _regionlist_to_single_region
from spectral_cube import SpectralCube
from spectral_cube.wcs_utils import strip_wcs_from_header
from spectral_cube.utils import NoBeamError
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy import log
from astropy.utils.console import ProgressBar
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
import warnings
import pathos
from spectral_cube.utils import SpectralCubeWarning
warnings.filterwarnings(action='ignore', category=SpectralCubeWarning,
                        append=True)


def read_as_2d(fn, minval=None):
    print(".", end='', flush=True)
    if fn.endswith('fits'):
        fh = fits.open(fn)
        ww = wcs.WCS(fh[0].header).celestial

        # strip the WCS
        header = strip_wcs_from_header(fh[0].header)
        header.update(ww.to_header())

        if minval is None:
            hdu2d = fits.PrimaryHDU(data=fh[0].data.squeeze(),
                                    header=header)
        else:
            data = fh[0].data.squeeze()
            # meant for weights; we're setting weights to zero
            data[data < minval] = 0
            # sanity check
            assert np.nanmax(data) == 1
            hdu2d = fits.PrimaryHDU(data=data,
                                    header=header)
        return fits.HDUList([hdu2d])
    else:
        cube = SpectralCube.read(fn, format='casa_image')
        return fits.HDUList([cube[0].hdu])


def get_peak(fn, slab_kwargs=None, rest_value=None):
    print(".", end='', flush=True)
    ft = 'fits' if fn.endswith(".fits") else "casa_image"
    cube = SpectralCube.read(fn, use_dask=True, format=ft).with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=rest_value)
    if slab_kwargs is not None:
        cube = cube.spectral_slab(**slab_kwargs)
    with cube.use_dask_scheduler('threads'):
        if cube.unit == u.dimensionless_unscaled:
            mx = cube.max(axis=0)
        else:
            mx = cube.max(axis=0).to(u.K)
    return mx


def get_m0(fn, slab_kwargs=None, rest_value=None):
    print(".", end='', flush=True)
    ft = 'fits' if fn.endswith(".fits") else "casa_image"
    cube = SpectralCube.read(fn, use_dask=True, format=ft).with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=rest_value)
    if slab_kwargs is not None:
        cube = cube.spectral_slab(**slab_kwargs)
    with cube.use_dask_scheduler('threads'):
        moment0 = cube.moment0(axis=0)
    moment0 = (moment0 * u.s / u.km).to(u.K,
                                        equivalencies=cube.beam.jtok_equiv(cube.with_spectral_unit(u.GHz).spectral_axis.mean())) * u.km / u.s
    return moment0


def make_mosaic(twod_hdus, name, norm_kwargs={}, slab_kwargs=None,
                weights=None,
                target_header=None,
                commonbeam=False,
                beams=None,
                rest_value=None, cbar_unit=None, array='7m', basepath='./'):

    if target_header is None:
        log.info(f"Finding WCS for {len(twod_hdus)} HDUs")
        target_wcs, shape_out = find_optimal_celestial_wcs(twod_hdus, frame='galactic')
    else:
        target_wcs = wcs.WCS(target_header)
        shape_out = (target_header['NAXIS2'], target_header['NAXIS1'])

    if commonbeam:
        if beams is None:
            beams = radio_beam.Beams(beams=[radio_beam.Beam.from_fits_header(hdul[0].header)
                                            for hdul in twod_hdus])
        if isinstance(commonbeam, radio_beam.Beam):
            cb = commonbeam
        else:
            cb = beams.common_beam()

        if isinstance(commonbeam, str) and commonbeam == 'circular':
            circbeam = radio_beam.Beam(major=cb.major, minor=cb.major, pa=0)
            cb = circbeam

        log.info("Loading HDUs and projecting to common beam")
        prjs = [spectral_cube.Projection.from_hdu(hdul) for hdul in
                ProgressBar(twod_hdus)]
        for prj, bm in (zip(ProgressBar(prjs), beams)):
            try:
                prj.beam
            except NoBeamError:
                prj.beam = bm

        log.info(f"Convolving HDUs to common beam {cb}\n")
        pb = ProgressBar(len(prjs))

        def reprj(prj):
            rslt = prj.convolve_to(cb).hdu
            pb.update()
            return rslt

        with pathos.Pool() as pool:
            twod_hdus = pool.map(reprj, prjs)

    log.info(f"Reprojecting and coadding {len(twod_hdus)} HDUs.\n")
    # number of items to count in progress bar
    npb = len(twod_hdus)
    pb = ProgressBar(npb)

    def repr_function(*args, **kwargs):
        rslt = reproject_interp(*args, **kwargs)
        pb.update()
        return rslt

    prjarr, footprint = reproject_and_coadd(twod_hdus,
                                            target_wcs,
                                            input_weights=weights,
                                            shape_out=shape_out,
                                            reproject_function=repr_function)
    header = target_wcs.to_header()
    if commonbeam:
        header.update(cb.to_header_keywords())

    log.info("Writing reprojected data to disk")
    fits.PrimaryHDU(data=prjarr, header=header).writeto(f'{basepath}/mosaics/{array}_{name}_mosaic.fits', overwrite=True)

    log.info("Creating plots")
    import pylab as pl
    from astropy import visualization
    pl.rc('axes', axisbelow=True)
    pl.matplotlib.use('agg')

    front = 10
    back = -10
    fronter = 20

    fig = pl.figure(figsize=(20, 7))
    ax = fig.add_subplot(111, projection=target_wcs)
    im = ax.imshow(prjarr, norm=visualization.simple_norm(prjarr, **norm_kwargs), zorder=front, cmap='gray')
    cbar = pl.colorbar(mappable=im)
    if cbar_unit is not None:
        cbar.set_label(cbar_unit)
    ax.coords[0].set_axislabel('Galactic Longitude')
    ax.coords[1].set_axislabel('Galactic Latitude')
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.coords[0].set_ticks(spacing=0.1 * u.deg)
    ax.coords[0].set_ticklabel(rotation=45, pad=20)

    fig.savefig(f'{basepath}/mosaics/{array}_{name}_mosaic.png', bbox_inches='tight')
    fig.savefig(f'{basepath}/mosaics/{array}_{name}_mosaic_hires.png', bbox_inches='tight', dpi=300)

    ax.coords.grid(True, color='black', ls='--', zorder=back)

    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='black', ls=':', zorder=back)
    overlay[0].set_axislabel('Right Ascension (ICRS)')
    overlay[1].set_axislabel('Declination (ICRS)')
    overlay[0].set_major_formatter('hh:mm')
    ax.set_axisbelow(True)
    ax.set_zorder(back)
    fig.savefig(f'{basepath}/mosaics/{array}_{name}_mosaic_withgrid.png', bbox_inches='tight')

    tbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/SB_naming.tsv', format='ascii.csv', delimiter='\t')

    log.info("Computing overlap regions")
    # create a list of composite regions for labeling purposes
    composites = []
    flagmap = np.zeros(prjarr.shape, dtype='int')
    # loop over the SB_naming table
    for row in ProgressBar(tbl):
        indx = row['Proposal ID'][3:]

        # load up regions
        regs = regions.Regions.read(f'{basepath}/reduction_ACES/aces/data/regions/final_cmz{indx}.reg')
        pregs = [reg.to_pixel(target_wcs) for reg in regs]

        composite = _regionlist_to_single_region(pregs)
        composite.meta['label'] = indx
        composites.append(composite)

        for comp in composites:
            cmsk = comp.to_mask()

            slcs_big, slcs_small = cmsk.get_overlap_slices(prjarr.shape)
            if slcs_big is not None and slcs_small is not None:
                flagmap[slcs_big] += (cmsk.data[slcs_small] * int(comp.meta['label'])) * (flagmap[slcs_big] == 0)
            else:
                print("Error - the composite mask has no overlap with the flag map.  "
                      "I don't know why this occurs but it definitely is not expected.")

    fits.PrimaryHDU(data=flagmap,
                    header=target_wcs.to_header()).writeto(f'{basepath}/mosaics/{array}_{name}_field_number_map.fits', overwrite=True)

    ax.contour(flagmap, cmap='prism', levels=np.arange(flagmap.max()) + 0.5, zorder=fronter)

    for ii in np.unique(flagmap):
        if ii > 0:
            fsum = (flagmap == ii).sum()
            cy, cx = ((np.arange(flagmap.shape[0])[:, None] * (flagmap == ii)).sum() / fsum,
                      (np.arange(flagmap.shape[1])[None, :] * (flagmap == ii)).sum() / fsum)
            pl.text(cx, cy, f"{ii}\n{tbl[ii-1]['Obs ID']}",
                    horizontalalignment='left', verticalalignment='center',
                    color=(1, 0.8, 0.5), transform=ax.get_transform('pixel'),
                    zorder=fronter)

    fig.savefig(f'{basepath}/mosaics/{array}_{name}_mosaic_withgridandlabels.png', bbox_inches='tight')

    if commonbeam:
        return cb
