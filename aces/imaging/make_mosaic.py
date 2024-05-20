import numpy as np
import regions
import radio_beam
from radio_beam.utils import BeamError
import spectral_cube
import PIL
from spectral_cube.lower_dimensional_structures import Projection
from spectral_cube.spectral_cube import _regionlist_to_single_region
from spectral_cube import SpectralCube
from spectral_cube.wcs_utils import strip_wcs_from_header
from spectral_cube.utils import NoBeamError
from spectral_cube.cube_utils import mosaic_cubes
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy.wcs import WCS
from astropy import log
from astropy.utils.console import ProgressBar
from astropy.convolution import convolve_fft, convolve, Gaussian2DKernel
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
import os
import glob
import copy
import shutil
from functools import partial
from multiprocessing import Process, Pool
import warnings
from spectral_cube.utils import SpectralCubeWarning

from tqdm.auto import tqdm

from aces import conf

basepath = conf.basepath

warnings.filterwarnings(action='ignore', category=SpectralCubeWarning,
                        append=True)


def get_common_beam(beams):
    for epsilon in (5e-4, 1e-3, 1e-4, 5e-3, 1e-2):
        for beam_threshold in np.logspace(-6, -2, 5):
            try:
                commonbeam = beams.common_beam(tolerance=beam_threshold, epsilon=epsilon)
                print(f"Successfully acquired common beam with tolerance={beam_threshold} and epsilon={epsilon}.  beam={commonbeam}")
                return commonbeam
            except (BeamError, ValueError) as ex:
                print(f"Encountered beam error '{ex}' with threshold {beam_threshold} and epsilon {epsilon}.  Trying again.")
    raise BeamError("Failed to find common beam.")


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
             folder=None, threshold=None, rel_threshold=None, fail_on_zeros=True):
    print(".", end='', flush=True)
    outfn = fn.replace(".fits", "") + f"{suffix}_max.fits"
    if folder is not None:
        outfn = os.path.join(folder, os.path.basename(outfn))
    if os.path.exists(outfn):
        hdu = fits.open(outfn)
        proj = Projection.from_hdu(hdu)
        if rel_threshold is not None:
            threshold = rel_threshold * proj.max()
            print(f"Set threshold to {threshold} based on rel_threshold={rel_threshold}")
        if threshold is not None:
            proj[proj < threshold] = 0
        if fail_on_zeros and np.nansum(proj) == 0:
            raise ValueError(f"File {fn} reduced to all zeros")
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
                    beam = get_common_beam(cube.beams)
                    equiv = beam.jtok_equiv(cube.with_spectral_unit(u.GHz).spectral_axis.mean())
                    mxjy = cube.max(axis=0)
                    if hasattr(mxjy, '_beam') and mxjy._beam is None:
                        mxjy._beam = beam
                    try:
                        assert hasattr(mxjy, 'beam')
                        assert mxjy.beam is not None
                    except Exception as ex:
                        print(ex)
                        mxjy = mxjy.with_beam(beam, raise_error_jybm=False)

                    mx = mxjy.to(u.K, equivalencies=equiv)
        if fail_on_zeros and np.nansum(mx.value) == 0:
            raise ValueError(f"File {fn} reduced to all zeros")
        if save_file:
            mx.hdu.writeto(outfn)

        if rel_threshold is not None:
            threshold = rel_threshold * mx.max()
            print(f"Set threshold to {threshold} based on rel_threshold={rel_threshold}")
        if threshold is not None:
            mx[mx.value < threshold] = 0
        return mx


def get_m0(fn, slab_kwargs=None, rest_value=None, suffix="", save_file=True, folder=None):
    print(".", end='', flush=True)
    outfn = fn.replace(".fits", "") + f"{suffix}_mom0.fits"
    if folder is not None:
        outfn = os.path.join(folder, os.path.basename(outfn))
    if os.path.exists(outfn):
        hdu = fits.open(outfn)
        proj = Projection.from_hdu(hdu)
        return proj
    else:
        ft = 'fits' if fn.endswith(".fits") else "casa_image"
        cube = SpectralCube.read(fn, use_dask=True, format=ft).with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=rest_value)
        cube.beam_threshold = 0.1  # SO2 or the one after it had 5% beam variance
        if slab_kwargs is not None:
            cube = cube.spectral_slab(**slab_kwargs)
        with cube.use_dask_scheduler('threads'):
            moment0 = cube.moment0(axis=0)
        if hasattr(cube, 'beam'):
            equiv = cube.beam.jtok_equiv(cube.with_spectral_unit(u.GHz).spectral_axis.mean())
        elif hasattr(cube, 'beams'):
            beam = get_common_beam(cube.beams)
            equiv = beam.jtok_equiv(cube.with_spectral_unit(u.GHz).spectral_axis.mean())

        moment0 = (moment0 * u.s / u.km).to(u.K,
                                            equivalencies=equiv) * u.km / u.s
        if save_file:
            moment0.hdu.writeto(outfn)
        return moment0


def check_hdus(hdus):
    bad = 0
    for hdu in hdus:
        if isinstance(hdu, fits.PrimaryHDU):
            ttl = np.nansum(hdu.data)
        else:
            ttl = np.nansum(hdu)
        if ttl == 0:
            bad = bad + 1

    if bad > 0:
        raise ValueError(f"Found {bad} bad HDUs (all zero)")


def makepng(data, wcs, imfn, footprint=None, **norm_kwargs):
    import pylab as pl
    from astropy import visualization
    import matplotlib.colors as mcolors
    import pyavm

    colors1 = pl.cm.gray_r(np.linspace(0., 1, 128))
    colors2 = pl.cm.hot(np.linspace(0, 1, 128))

    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    sel = np.isnan(data) | (data == 0) | (footprint == 0 if footprint is not None else False)

    norm = visualization.simple_norm(data, **norm_kwargs)
    if 'min_cut' in norm_kwargs:
        norm.vmin = norm_kwargs['min_cut']
        assert norm.vmin == norm_kwargs['min_cut']

    colordata = mymap(norm(data))
    ct = (colordata[:, :, :] * 256).astype('uint8')
    ct[(colordata[:, :, :] * 256) > 255] = 255
    ct[:, :, 3][sel] = 0
    img = PIL.Image.fromarray(ct[::-1, :, :])
    img.save(imfn)
    avm = pyavm.AVM.from_wcs(wcs, shape=data.shape)
    avm.embed(imfn, imfn)


def make_mosaic(twod_hdus, name, norm_kwargs={}, slab_kwargs=None,
                weights=None,
                target_header=None,
                commonbeam=None,
                beams=None,
                rest_value=None,
                cbar_unit=None,
                array='7m',
                folder=None,  # must be specified though
                basepath='./'):
    """
    Given a long list of 2D HDUs and an output name, make a giant mosaic.
    """

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
            cb = get_common_beam(beams)

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
    if commonbeam is not None:
        header.update(cb.to_header_keywords())

    assert not np.all(np.isnan(prjarr))

    outfile = f'{basepath}/mosaics/{folder}/{array}_{name}_mosaic.fits'
    log.info(f"Writing reprojected data to {outfile}")
    fits.PrimaryHDU(data=prjarr, header=header).writeto(outfile, overwrite=True)

    log.info("Creating plots")
    import pylab as pl
    from astropy import visualization
    import matplotlib.colors as mcolors
    import pyavm

    colors1 = pl.cm.gray_r(np.linspace(0., 1, 128))
    colors2 = pl.cm.hot(np.linspace(0, 1, 128))

    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

    pl.rc('axes', axisbelow=True)
    pl.matplotlib.use('agg')

    front = 10
    back = -10
    fronter = 20

    fig = pl.figure(figsize=(20, 7))
    ax = fig.add_subplot(111, projection=target_wcs)
    norm = visualization.simple_norm(prjarr, **norm_kwargs)
    im = ax.imshow(prjarr, norm=norm, zorder=front, cmap=mymap)
    cbar = pl.colorbar(mappable=im)
    if cbar_unit is not None:
        cbar.set_label(cbar_unit)
    ax.coords[0].set_axislabel('Galactic Longitude')
    ax.coords[1].set_axislabel('Galactic Latitude')
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.coords[0].set_ticks(spacing=0.1 * u.deg)
    ax.coords[0].set_ticklabel(rotation=45, pad=20)

    fig.savefig(f'{basepath}/mosaics/{folder}/{array}_{name}_mosaic.png', bbox_inches='tight')
    fig.savefig(f'{basepath}/mosaics/{folder}/{array}_{name}_mosaic_hires.png', bbox_inches='tight', dpi=300)

    ax.coords.grid(True, color='black', ls='--', zorder=back)

    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='black', ls=':', zorder=back)
    overlay[0].set_axislabel('Right Ascension (ICRS)')
    overlay[1].set_axislabel('Declination (ICRS)')
    overlay[0].set_major_formatter('hh:mm')
    ax.set_axisbelow(True)
    ax.set_zorder(back)
    fig.savefig(f'{basepath}/mosaics/{folder}/{array}_{name}_mosaic_withgrid.png', bbox_inches='tight')

    tbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/SB_naming.md', format='ascii.csv', delimiter='|', data_start=2)

    log.info("Creating AVM-embedded colormapped image")

    imfn = f'{basepath}/mosaics/{folder}/{array}_{name}_mosaic_noaxes.png'
    makepng(prjarr, target_wcs, imfn, footprint=footprint, **norm_kwargs)

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

    outfile = f'{basepath}/mosaics/{folder}/{array}_{name}_field_number_map.fits'
    log.info(f"Writing flag image to {outfile}")
    fits.PrimaryHDU(data=flagmap,
                    header=target_wcs.to_header()).writeto(outfile, overwrite=True)

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

    fig.savefig(f'{basepath}/mosaics/{folder}/{array}_{name}_mosaic_withgridandlabels.png', bbox_inches='tight')

    # close figures to avoid memory leak
    pl.close('all')

    if commonbeam is not None:
        return cb


def all_lines(header, parallel=False, array='12m', glob_suffixes=('cube.I.iter1.image.pbcor', 'cube.I.manual.image.pbcor'),
              lines='all', folder='',
              globdir='calibrated/working', use_weights=True):

    from astropy.table import Table

    tbl = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/linelist.csv')

    if parallel:
        processes = []

    for row in tbl:
        spwn = row[f'{array} SPW']
        line = row['Line'].replace(" ", "_").replace("(", "_").replace(")", "_")
        restf = row['Rest (GHz)'] * u.GHz

        if lines != 'all' and not (line in lines or row['Line'] in lines):
            continue

        log.info(f"{array} {line} {restf}")

        filelist = []
        for glob_suffix in glob_suffixes:
            filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/{globdir}/*spw{spwn}.{glob_suffix}')
        assert len(filelist) > 0
        if use_weights:
            weightfiles = [fn.replace("image.pbcor", "weight") for fn in filelist]
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
            check_hdus(hdus)
            if use_weights:
                wthdus = pool.map(partial(get_peak,
                                          **{'slab_kwargs': {'lo': -2 * u.km / u.s, 'hi': 2 * u.km / u.s},
                                             'rest_value': restf},
                                          suffix=f'_{line}',
                                          rel_threshold=0.25,  # pb limit
                                          ),
                                  weightfiles)
                wthdus = [x.hdu for x in wthdus]
                check_hdus(wthdus)
        else:
            hdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s},
                             rest_value=restf, suffix=f'_{line}').hdu for fn in filelist]
            check_hdus(hdus)
            print(flush=True)
            if use_weights:
                wthdus = [get_peak(fn, slab_kwargs={'lo': -2 * u.km / u.s, 'hi': 2 * u.km / u.s},
                                   rest_value=restf, suffix=f'_{line}',
                                   rel_threshold=0.25,  # pb limit
                                   ).hdu for fn in weightfiles]
                check_hdus(wthdus)
                print(flush=True)

        if parallel:
            print(f"Starting function make_mosaic for {line} peak intensity (parallel mode)")
            proc = Process(target=make_mosaic, args=(hdus, ),
                           kwargs=dict(name=f'{line}_max', cbar_unit='K',
                           norm_kwargs=dict(max_cut=5, min_cut=-0.01, stretch='asinh'),
                           array=array, basepath=basepath, weights=wthdus if use_weights else None,
                           folder=folder,
                           target_header=header))
            proc.start()
            processes.append(proc)
        else:
            print(f"Starting function make_mosaic for {line} peak intensity")
            make_mosaic(hdus, name=f'{line}_max', cbar_unit='K',
                        norm_kwargs=dict(max_cut=5, min_cut=-0.01, stretch='asinh'),
                        array=array, basepath=basepath, weights=wthdus if use_weights else None,
                        folder=folder,
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
            print(f"Starting function make_mosaic for {line} moment 0")
            proc = Process(target=make_mosaic, args=(m0hdus, ),
                           kwargs=dict(name=f'{line}_m0', cbar_unit='K km/s',
                           norm_kwargs=dict(max_cut=20, min_cut=-1, stretch='asinh'),
                           folder=folder,
                           array=array, basepath=basepath, weights=wthdus if use_weights else None, target_header=header))
            proc.start()
            processes.append(proc)
        else:
            make_mosaic(m0hdus, name=f'{line}_m0', cbar_unit='K km/s',
                        norm_kwargs=dict(max_cut=20, min_cut=-1, stretch='asinh'),
                        array=array, basepath=basepath,
                        folder=folder,
                        weights=wthdus if use_weights else None, target_header=header)
            del m0hdus
            if use_weights:
                del wthdus

    if parallel:
        for proc in processes:
            proc.join()


def make_giant_mosaic_cube_header(target_header,
                                  reference_frequency,
                                  cdelt_kms,
                                  nchan,
                                  test=False,
                                  ):
    header = fits.Header.fromtextfile(target_header)
    header['NAXIS'] = 3
    header['CDELT3'] = cdelt_kms
    header['CUNIT3'] = 'km/s'
    header['CRVAL3'] = 0
    header['NAXIS3'] = 3 if test else nchan
    header['CRPIX3'] = header['NAXIS3'] // 2
    header['CTYPE3'] = 'VRAD'
    restfrq = u.Quantity(reference_frequency, u.Hz).to(u.Hz).value
    header['RESTFRQ'] = restfrq
    header['SPECSYS'] = 'LSRK'

    return header


def make_giant_mosaic_cube_channels(header, cubes, weightcubes, commonbeam,
                                    cubename,
                                    verbose=True,
                                    working_directory='/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/',
                                    channelmosaic_directory=f'{basepath}/mosaics/HNCO_Channels/',
                                    fail_if_cube_dropped=True,
                                    channels='all'):
    ww = WCS(header)
    wws = ww.spectral

    if not os.path.exists(channelmosaic_directory):
        os.mkdir(channelmosaic_directory)

    # Part 5: Make per-channel mosaics
    if verbose:
        print(f"Making per-channel mosaics into cube {cubename} with channels {channels}")
    pbar = tqdm(channels, desc='Channels (mosaic)') if verbose else channels
    for chan in pbar:
        # For each channel, we're making a full-frame mosaic
        theader = header.copy()

        theader['NAXIS3'] = 1
        theader['CRPIX3'] = 1
        theader['CRVAL3'] = wws.pixel_to_world(chan).to(u.km / u.s).value

        if verbose:
            pbar.set_description(f'Channel {chan}')

        chanfn = f'{working_directory}/{cubename}_CubeMosaic_channel{chan}.fits'
        if os.path.exists(f'{channelmosaic_directory}/{os.path.basename(chanfn)}'):
            if check_channel(f'{channelmosaic_directory}/{os.path.basename(chanfn)}', verbose=verbose):
                print(f"Skipping completed channel {chan}", flush=True)
            else:
                print(f"Removing failed channel {chan}.  Will need to re-run", flush=True)
                os.remove(f'{channelmosaic_directory}/{os.path.basename(chanfn)}')
                raise ValueError("Empty channel found - re-run needed.")
        elif os.path.exists(chanfn):
            if check_channel(chanfn, verbose=True):
                print(f"Channel {chan} was already completed, moving it to {channelmosaic_directory}", flush=True)
                shutil.move(chanfn, channelmosaic_directory)
            else:
                print(f"Removing failed channel {chan}={chanfn}.  Will need to re-run", flush=True)
                os.remove(chanfn)
                raise ValueError("Empty channel found - re-run needed.")
        else:
            print(f"Starting mosaic_cubes for channel {chan}", flush=True)
            #cubes=cubes[:2] # TEMP DEBUG
            #weightcubes=weightcubes[:2] # TEMP DEBUG
            mosaic_cubes(cubes,
                         target_header=theader,
                         commonbeam=commonbeam,
                         weightcubes=weightcubes,
                         spectral_block_size=2,
                         output_file=chanfn,
                         method='channel',
                         verbose=verbose,
                         fail_if_cube_dropped=fail_if_cube_dropped,
                         )
            print(f"\nChannel {chan} appears to have completed successfully, but we're checking first.", flush=True)

            if not np.any(np.isnan(fits.getdata(chanfn))):
                print(f"Channel {chan} had no nans but {(fits.getdata(chanfn) == 0).sum()} zeros.  Setting 0->nan")
                print("This is a bug that appeared sometime in Jan 2024 and I haven't been able to ID a cause")
                fh = fits.open(chanfn, mode='update')
                fh[0].data[fh[0].data == 0] = np.nan
                fh.close()

            if not check_channel(chanfn, verbose=verbose):
                raise ValueError("Produced an empty channel - raising exception as this is not expected")
            else:
                print(f"Channel {chan} completed successfully, moving it to {channelmosaic_directory}", flush=True)
                shutil.move(chanfn, channelmosaic_directory)

        if verbose:
            print("\n\n", flush=True)


def check_channel(chanfn, verbose=True):
    data = fits.getdata(chanfn)
    if np.all(np.isnan(data)) or np.nansum(data) == 0:
        if verbose:
            print(f"{chanfn} failed: data.sum={data.sum()}, nansum={np.nansum(data)} data.std={data.std()} data.finite={np.sum(~np.isnan(data))}")
        return False
    elif (not np.any(np.isnan(data))) and (np.nansum(data == 0) > 1000):
        if verbose:
            print(f"{chanfn} failed: data.sum={data.sum()}, nansum={np.nansum(data)} data.std={data.std()} data.finite={np.sum(~np.isnan(data))} data.notfinite={np.sum(np.isnan(data))}")
        return False
    else:
        if verbose:
            print(f"{chanfn} succeeded: data.sum={data.sum()}, nansum={np.nansum(data)} data.std={data.std()} data.finite={np.sum(~np.isnan(data))}")
        return True


def make_giant_mosaic_cube(filelist,
                           reference_frequency,
                           cdelt_kms,
                           cubename,
                           nchan,
                           test=False, verbose=True,
                           weightfilelist=None,
                           target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr',
                           working_directory='/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/',
                           channelmosaic_directory=f'{basepath}/mosaics/HNCO_Channels/',
                           beam_threshold=3.2 * u.arcsec,
                           channels='all',
                           fail_if_cube_dropped=True,
                           skip_channel_mosaicing=False,
                           skip_final_combination=False,
                           ):
    """
    This takes too long as a full cube, so we have to do it slice-by-slice

    channels : 'all' or a list of ints
        This gives you the option to run only one channel at a time
    beam_threshold : angle-like quantity
        Cubes with beams larger than this will be excldued
    """

    reference_frequency = u.Quantity(reference_frequency, u.Hz)

    # Part 1: Make the header
    header = make_giant_mosaic_cube_header(target_header=target_header,
                                           reference_frequency=reference_frequency,
                                           cdelt_kms=cdelt_kms,
                                           nchan=nchan,
                                           test=test)

    # Part 2: Load the cubes
    print("Converting spectral units", flush=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cubes = [SpectralCube.read(fn,
                                   format='fits' if fn.endswith('fits') else 'casa_image',
                                   use_dask=True).with_spectral_unit(u.km / u.s,
                                                                     velocity_convention='radio',
                                                                     rest_value=reference_frequency)
                 for fn in filelist]
        if weightfilelist is None:
            weightcubes = [(SpectralCube.read(fn.replace(".image.pbcor", ".weight"),
                                              format='fits' if fn.endswith('fits') else 'casa_image', use_dask=True)
                            .with_spectral_unit(u.km / u.s,
                                                velocity_convention='radio',
                                                rest_value=reference_frequency)
                            ) for fn in filelist]
        else:
            weightcubes = [SpectralCube.read(fn, format='fits' if fn.endswith('fits') else 'casa_image',
                                             use_dask=True)
                           .with_spectral_unit(u.km / u.s,
                                               velocity_convention='radio',
                                               rest_value=reference_frequency)
                           for fn in weightfilelist]

    # BUGFIX: there are FITS headers that incorrectly specify UTC in caps
    for cube in cubes + weightcubes:
        cube._wcs.wcs.timesys = cube.wcs.wcs.timesys.lower()
        if hasattr(cube.mask, '_wcs'):
            cube.mask._wcs.wcs.timesys = cube.wcs.wcs.timesys.lower()

    # Part 3: Filter out bad cubes
    # flag out wild outliers
    # there are 2 as of writing
    if verbose:
        print("Filtering out cubes with sketchy beams", flush=True)
    for cube, fn in zip(cubes, filelist):
        try:
            cube.beam if hasattr(cube, 'beam') else cube.beams
        except NoBeamError as ex:
            print(f"{fn} has no beam")
            raise ex
    beams = [get_common_beam(cube.beams) if hasattr(cube, 'beams') else cube.beam
             for cube in cubes]
    ok = [beam.major < beam_threshold for beam in beams]
    if verbose:
        if not all(ok):
            print(f"Filtered down to {np.sum(ok)} of {len(ok)} cubes with beam major > {beam_threshold}")
            print(f"Filtered cubes include: {[fn for fn, k in zip(filelist, ok) if not k]}")
        else:
            print(f"Found {np.sum(ok)} cubes with good beams (i.e., all of them)")

    cubes = [cube for k, cube in zip(ok, cubes) if k]
    weightcubes = [cube for k, cube in zip(ok, weightcubes) if k]

    if verbose:
        print(f"There are {len(cubes)} cubes and {len(weightcubes)} weightcubes.", flush=True)

    # Part 4: Determine common beam
    if verbose:
        print("Determining common beam")
    beams = radio_beam.Beams(beams=[cube.beam
                                    if hasattr(cube, 'beam')
                                    else get_common_beam(cube.beams)
                                    for cube in cubes])
    commonbeam = get_common_beam(beams)
    header.update(commonbeam.to_header_keywords())

    if channels == 'all':
        channels = range(header['NAXIS3'])
    elif channels == 'slurm':
        channels = slurm_set_channels(nchan)

    # Part 5: make per-channel mosaics
    if not skip_channel_mosaicing:
        if verbose:
            print(f"Making the channels from our {len(cubes)} cubes and {len(weightcubes)} weightcubes for channels {channels}")
        make_giant_mosaic_cube_channels(header, cubes, weightcubes,
                                        commonbeam=commonbeam,
                                        cubename=cubename,
                                        verbose=verbose,
                                        working_directory=working_directory,
                                        channelmosaic_directory=channelmosaic_directory,
                                        fail_if_cube_dropped=fail_if_cube_dropped,
                                        channels=channels,
                                        )
    else:
        print("Skipped channel mosaicking")

    if not skip_final_combination and not test:
        if verbose:
            print("Combining channels into mosaic cube")
        combine_channels_into_mosaic_cube(header,
                                          cubename,
                                          nchan,
                                          channels=channels,
                                          working_directory=working_directory,
                                          channelmosaic_directory=channelmosaic_directory,
                                          verbose=verbose,
                                          )


def combine_channels_into_mosaic_cube(header, cubename, nchan, channels,
                                      working_directory='/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/',
                                      channelmosaic_directory=f'{basepath}/mosaics/HNCO_Channels/',
                                      verbose=False,
                                      ):
    # Part 6: Create output supergiant cube into which final product will be stashed
    output_working_file = f'{working_directory}/{cubename}_CubeMosaic.fits'
    output_file = f'{basepath}/mosaics/cubes/{cubename}_CubeMosaic.fits'
    if verbose:
        print(f"Beginning combination: working file is {output_working_file} (exists={os.path.exists(output_working_file)}) and final output is {output_file} (exists={os.path.exists(output_file)})")
    if os.path.exists(output_working_file) and not os.path.exists(output_file):
        if verbose:
            print(f"Continuing work on existing file {output_working_file}")
    elif not os.path.exists(output_file):
        if verbose:
            print(f"Making new file {output_file}")
        # Make a new file
        # https://docs.astropy.org/en/stable/generated/examples/io/skip_create-large-fits.html#sphx-glr-generated-examples-io-skip-create-large-fits-py
        hdu = fits.PrimaryHDU(data=np.ones([5, 5, 5], dtype=float),
                              header=header,
                              )
        for kwd in ('NAXIS1', 'NAXIS2', 'NAXIS3'):
            hdu.header[kwd] = header[kwd]

        assert header['NAXIS3'] == nchan, f"Sanity check: header must have same number of channels as requested (header={header['NAXIS3']}, nchan={nchan})"
        shape_opt = header['NAXIS3'], header['NAXIS2'], header['NAXIS1']

        if verbose:
            print(f"Creating output working file {output_working_file}")
        header.tofile(output_working_file, overwrite=True)
        with open(output_working_file, 'rb+') as fobj:
            fobj.seek(len(header.tostring()) +
                      (np.prod(shape_opt) * np.abs(header['BITPIX'] // 8)) - 1)
            fobj.write(b'\0')
    elif os.path.exists(output_file) and not os.path.exists(output_working_file):
        if verbose:
            print(f"Working on file {output_file}, but moving it to {output_working_file} first")
        shutil.move(output_file, output_working_file)
        if verbose:
            print(f"Completed move of {output_file} to {output_working_file}")
    else:
        raise ValueError(f"This outcome should not be possible: both {output_file} and {output_working_file} exist")

    hdu = fits.open(output_working_file)
    if hdu[0].data.shape[0] != nchan:
        raise ValueError(f"Existing file on disk {output_working_file} has {hdu[0].data.shape[0]} channels instead of the requested {nchan}")

    hdu = fits.open(output_working_file, mode='update', overwrite=True)
    output_array = hdu[0].data
    hdu.flush()  # make sure the header gets written right

    # Part 7: Populate supergiant cube by copying data over, channel-by-channel
    if verbose:
        print(f"Beginning channel filling from channel mosaic directory {channelmosaic_directory}")
        pbar = tqdm(channels, desc='Channels')
    else:
        pbar = channels
    for chan in pbar:
        chanfn = f'{channelmosaic_directory}/{cubename}_CubeMosaic_channel{chan}.fits'
        if os.path.exists(chanfn):
            if verbose:
                pbar.set_description('Channels (filling)')
            output_array[chan] = fits.getdata(chanfn)
            if verbose:
                pbar.set_description('Channels (flushing)')
            hdu.flush()

    if verbose:
        print(f"Moving {output_working_file} to {output_file}")
    shutil.move(output_working_file, output_file)
    assert os.path.exists(output_file), f"Failed to move {output_working_file} to {output_file}"
    if verbose:
        print(f"Successfully moved {output_working_file} to {output_file}")


def slurm_set_channels(nchan):
    if os.getenv('SLURM_ARRAY_TASK_ID') is not None:
        slurm_array_task_id = int(os.getenv('SLURM_ARRAY_TASK_ID'))
        slurm_array_task_count = int(os.getenv('SLURM_ARRAY_TASK_COUNT'))

        nchan_per = nchan / slurm_array_task_count
        if nchan_per != int(nchan_per):
            raise ValueError(f"Must divide up slurm jobs evenly.  Got {nchan_per} channels, which isn't an int")
        nchan_per = int(nchan_per)
        channels = list(range(slurm_array_task_id * nchan_per,
                              (slurm_array_task_id + 1) * nchan_per))
        return channels


def make_downsampled_cube(cubename, outcubename, factor=9, overwrite=True, use_dask=True, spectrally_too=True):
    """
    TODO: may need to dump-to-temp while reprojecting
    """
    cube = SpectralCube.read(cubename, use_dask=use_dask)
    if use_dask:
        import dask.array as da
        from dask.diagnostics import ProgressBar
    else:
        cube.allow_huge_operations = True
        import contextlib
        ProgressBar = contextlib.nullcontext

    print(f"Downsampling cube {cubename} -> {outcubename}")
    print(cube)
    start = 0
    with ProgressBar():
        dscube = cube[:, start::factor, start::factor]
        # this is a hack; see https://github.com/astropy/astropy/pull/10897
        # the spectral-cube approach is right normally, but we're cheating here and just taking every 9th pixel, we're not smoothing first.
        dscube.wcs.celestial.wcs.crpix[:] = (dscube.wcs.celestial.wcs.crpix[:] - 1 - start) / factor + 1
        dscube.write(outcubename, overwrite=overwrite)
    # else:
    #     cube = SpectralCube.read(cubename, use_dask=False)
    #     #hdr = cube.header.copy()
    #     #hdr['NAXIS1'] = int(hdr['NAXIS1'] / factor)
    #     #hdr['NAXIS2'] = int(hdr['NAXIS2'] / factor)
    #     #hdr['CRPIX1'] = int(hdr['CRPIX1'] / factor)
    #     #hdr['CRPIX2'] = int(hdr['CRPIX2'] / factor)
    #     #hdr['CDELT1'] = (hdr['CDELT1'] * factor)
    #     #hdr['CDELT2'] = (hdr['CDELT2'] * factor)
    #     cube.allow_huge_operations = True
    #     print(f"Downsampling cube {cubename} -> {outcubename}")
    #     print(cube)
    #     dscube = cube[:, ::factor, ::factor]
    #     dscube.write(outcubename, overwrite=overwrite)
    #         #dscube = cube.reproject(hdr, return_type='dask', filled=False, parallel=True, block_size=[1,1000,1000])
    #         ## non-dask dscube = cube.reproject(hdr, use_memmap=True)
    #         #print("Writing")
    #         #dscube.write(outcubename, overwrite=overwrite)

    if spectrally_too:
        dscube_s = dscube.downsample_axis(factor=factor, axis=0)
        assert outcubename.endswith('.fits')
        dscube_s.write(outcubename.replace(".fits", "_spectrally.fits"), overwrite=True)


def rms_map(img, kernel=Gaussian2DKernel(10)):
    """
    Gaussian2DKernel should be larger than the beam, probably at least 2x larger
    """
    # nans = np.isnan(img)
    sm = convolve_fft(img, kernel, allow_huge=True)
    res = img - sm
    var = res**2
    smvar = convolve_fft(var, kernel, allow_huge=True)
    rms = smvar**0.5
    # restore NaNs: the convolution process will naturally fill them
    # don't do this: it turns internal pixels into nans that shouldn't be
    # rms[nans] = np.nan
    return rms


def rms(prefix='12m_continuum', folder='12m_flattened', threshold=2.5, nbeams=3, maxiter=50):
    for fn in glob.glob(f'{basepath}/mosaics/{folder}/{prefix}*mosaic.fits'):
        if 'rms' in fn:
            # don't re-rms rms maps
            continue
        print(f"Computing RMS map for {fn}", flush=True)
        fh = fits.open(fn)
        ww = WCS(fh[0].header)
        pixscale = ww.proj_plane_pixel_area()**0.5
        try:
            beam = radio_beam.Beam.from_fits_header(fh[0].header)
            kernelwidth = ((beam.major * nbeams) / pixscale).decompose()
        except Exception as ex:
            print(ex, fn)
            kernelwidth = (2.5 * u.arcsec / pixscale).decompose()

        nans = np.isnan(fh[0].data)
        rms = rms_map(fh[0].data, kernel=Gaussian2DKernel(kernelwidth))
        rms[nans] = np.nan

        outname = fn.replace("_mosaic.fits", "_rms_mosaic.fits")
        fits.PrimaryHDU(data=rms, header=fh[0].header).writeto(outname, overwrite=True)
        print(f"Finished RMS map for {fn}")

        datacopy = copy.copy(fh[0].data)
        ndet = []
        for ii in range(maxiter):
            # second iteration: threshold and try again
            detections = (datacopy / rms) > threshold
            ndet_this = detections.sum()
            if ndet_this == 0:
                print(f"Converged in {ii} iterations")
                break
            ndet.append(ndet_this)
            print(f"Iteration {ii} detected {ndet}")

            datacopy[detections] = np.nan
            rms = rms_map(datacopy, kernel=Gaussian2DKernel(kernelwidth))
        rms[nans] = np.nan

        outname = fn.replace("_mosaic.fits", "_maskedrms_mosaic.fits")
        fits.PrimaryHDU(data=rms, header=fh[0].header).writeto(outname, overwrite=True)
        print(f"Finished masked RMS map for {fn}.  Detections iterated: {ndet}")
