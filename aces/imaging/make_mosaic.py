import numpy as np
import regions
import radio_beam
import spectral_cube
from spectral_cube.lower_dimensional_structures import Projection
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
import os
import glob
from functools import partial
from multiprocessing import Process, Pool
import warnings
from spectral_cube.utils import SpectralCubeWarning

from aces import conf

basepath = conf.basepath

warnings.filterwarnings(action='ignore', category=SpectralCubeWarning,
                        append=True)


def read_as_2d(fn, minval=None):
    print(".", end='', flush=True)
    if fn.endswith('fits') or fn.endswith('.fits.gz'):
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


def get_peak(fn, slab_kwargs=None, rest_value=None, suffix="", save_file=True,
             threshold=None):
    print(".", end='', flush=True)
    outfn = fn.replace(".fits", "") + f"{suffix}_max.fits"
    if os.path.exists(outfn):
        hdu = fits.open(outfn)
        proj = Projection.from_hdu(hdu)
        if threshold is not None:
            proj[proj < threshold] = 0
        return proj
    else:
        ft = 'fits' if fn.endswith(".fits") else "casa_image"
        cube = SpectralCube.read(fn, use_dask=True, format=ft).with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=rest_value)
        if slab_kwargs is not None:
            cube = cube.spectral_slab(**slab_kwargs)
        with cube.use_dask_scheduler('threads'):
            if cube.unit == u.dimensionless_unscaled:
                mx = cube.max(axis=0)
            else:
                if hasattr(cube, 'beam'):
                    mx = cube.max(axis=0).to(u.K)
                else:
                    log.warn(f"File {fn} is a multi-beam cube.")
                    beam = cube.beams.common_beam()
                    mx = cube.max(axis=0).to(u.K, beam.jtok_equiv(cube.spectral_axis.mean()))
        if save_file:
            mx.hdu.writeto(outfn)
        if threshold is not None:
            mx[mx.value < threshold] = 0
        return mx


def get_m0(fn, slab_kwargs=None, rest_value=None, suffix="", save_file=True):
    print(".", end='', flush=True)
    outfn = fn.replace(".fits", "") + f"{suffix}_mom0.fits"
    if os.path.exists(outfn):
        hdu = fits.open(outfn)
        proj = Projection.from_hdu(hdu)
        return proj
    else:
        ft = 'fits' if fn.endswith(".fits") else "casa_image"
        cube = SpectralCube.read(fn, use_dask=True, format=ft).with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=rest_value)
        if slab_kwargs is not None:
            cube = cube.spectral_slab(**slab_kwargs)
        with cube.use_dask_scheduler('threads'):
            moment0 = cube.moment0(axis=0)
        moment0 = (moment0 * u.s / u.km).to(u.K,
                                            equivalencies=cube.beam.jtok_equiv(cube.with_spectral_unit(u.GHz).spectral_axis.mean())) * u.km / u.s
        if save_file:
            moment0.hdu.writeto(outfn)
        return moment0


def make_mosaic(twod_hdus, name, norm_kwargs={}, slab_kwargs=None,
                weights=None,
                target_header=None,
                commonbeam=None,
                beams=None,
                rest_value=None,
                cbar_unit=None,
                array='7m',
                basepath='./'):

    if target_header is None:
        log.info(f"Finding WCS for {len(twod_hdus)} HDUs")
        target_wcs, shape_out = find_optimal_celestial_wcs(twod_hdus, frame='galactic')
    else:
        target_wcs = wcs.WCS(target_header)
        shape_out = (target_header['NAXIS2'], target_header['NAXIS1'])

    if commonbeam is not None:
        if beams is None:
            beams = radio_beam.Beams(beams=[radio_beam.Beam.from_fits_header(hdul[0].header)
                                            for hdul in twod_hdus])
            if array == '12m':
                for beam in beams:
                    assert beam.major < 3 * u.arcsec
        if isinstance(commonbeam, radio_beam.Beam):
            cb = commonbeam
        else:
            cb = beams.common_beam()

        if isinstance(commonbeam, str) and commonbeam == 'circular':
            circbeam = radio_beam.Beam(major=cb.major, minor=cb.major, pa=0)
            cb = circbeam
            if array == '12m':
                assert cb.major < 3 * u.arcsec

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

        # parallelization of this failed
        twod_hdus = [reprj(prj) for prj in prjs]

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

    outfile = f'{basepath}/mosaics/{array}_{name}_mosaic.fits'
    log.info(f"Writing reprojected data to {outfile}")
    fits.PrimaryHDU(data=prjarr, header=header).writeto(outfile, overwrite=True)

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

    tbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/SB_naming.md', format='ascii.csv', delimiter='|', data_start=2)

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

    outfile = f'{basepath}/mosaics/{array}_{name}_field_number_map.fits'
    log.info(f"Writing flag image to {outfile}")
    fits.PrimaryHDU(data=flagmap,
                    header=target_wcs.to_header()).writeto(outfile, overwrite=True)

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


def all_lines(header, parallel=False, array='12m', glob_suffix='cube.I.iter1.image.pbcor',
              globdir='calibrated/working', use_weights=True):

    from astropy.table import Table

    tbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/linelist.csv')

    if parallel:
        processes = []

    for row in tbl:
        spwn = row[f'{array} SPW']
        line = row['Line'].replace(" ", "_").replace("(", "_").replace(")", "_")
        restf = row['Rest (GHz)'] * u.GHz

        log.info(f"{array} {line}")

        filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/{globdir}/*spw{spwn}.{glob_suffix}')
        assert len(filelist) > 0
        if use_weights:
            weightfiles = [fn.replace("image.pbcor", "pb") for fn in filelist]
            for ii, (ifn, wfn) in enumerate(zip(list(filelist), list(weightfiles))):
                if not os.path.exists(wfn):
                    log.error(f"Missing file {wfn}")
                    filelist.remove(ifn)
                    weightfiles.remove(wfn)

        if parallel:
            pool = Pool()
            hdus = pool.map(partial(get_peak,
                                    **{'slab_kwargs': {'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s},
                                       'rest_value': restf},
                                    suffix=f'_{line}',
                                    ),
                            filelist,
                            )
            hdus = [x.hdu for x in hdus]
            if use_weights:
                wthdus = pool.map(partial(get_peak,
                                          **{'slab_kwargs': {'lo': -2 * u.km / u.s, 'hi': 2 * u.km / u.s},
                                             'rest_value': restf},
                                          suffix=f'_{line}',
                                          threshold=0.5,  # pb limit
                                          ),
                                  weightfiles)
                wthdus = [x.hdu for x in wthdus]
        else:
            hdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s},
                             rest_value=restf, suffix=f'_{line}').hdu for fn in filelist]
            print(flush=True)
            if use_weights:
                wthdus = [get_peak(fn, slab_kwargs={'lo': -2 * u.km / u.s, 'hi': 2 * u.km / u.s},
                                   rest_value=restf, suffix=f'_{line}',
                                   threshold=0.5,  # pb limit
                                   ).hdu for fn in weightfiles]
                print(flush=True)

        if parallel:
            print(f"Starting function make_mosaic for {line}")
            proc = Process(target=make_mosaic, args=(hdus,),
                           kwargs=dict(name=f'{line}_max', cbar_unit='K',
                           norm_kwargs=dict(max_cut=5, min_cut=-0.01, stretch='asinh'),
                           array=array, basepath=basepath, weights=wthdus if use_weights else None,
                           target_header=header))
            proc.start()
            processes.append(proc)
        else:
            make_mosaic(hdus, name=f'{line}_max', cbar_unit='K',
                        norm_kwargs=dict(max_cut=5, min_cut=-0.01, stretch='asinh'),
                        array=array, basepath=basepath, weights=wthdus if use_weights else None,
                        target_header=header)
            del hdus

        if parallel:
            pool = Pool()
            m0hdus = pool.map(partial(get_m0, **{'slab_kwargs': {'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, 'rest_value': restf}, suffix=f'_{line}'), filelist)
            m0hdus = [x.hdu for x in m0hdus]
        else:
            m0hdus = [get_m0(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=restf, suffix=f'_{line}').hdu for fn in filelist]
            print(flush=True)

        if parallel:
            print(f"Starting function make_mosaic for {line}")
            proc = Process(target=make_mosaic, args=(m0hdus,),
                           kwargs=dict(name=f'{line}_m0', cbar_unit='K km/s',
                           norm_kwargs=dict(max_cut=20, min_cut=-1, stretch='asinh'),
                           array=array, basepath=basepath, weights=wthdus if use_weights else None, target_header=header))
            proc.start()
            processes.append(proc)
        else:
            make_mosaic(m0hdus, name=f'{line}_m0', cbar_unit='K km/s',
                        norm_kwargs=dict(max_cut=20, min_cut=-1, stretch='asinh'),
                        array=array, basepath=basepath, weights=wthdus if use_weights else None, target_header=header)
            del m0hdus
            if use_weights:
                del wthdus

    if parallel:
        for proc in processes:
            proc.join()
