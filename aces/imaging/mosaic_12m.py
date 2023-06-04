"""
After running this, symlink to web via:
ln -s /orange/adamginsburg/ACES/mosaics/*m0_mosaic*png /orange/adamginsburg/web/secure/ACES/mosaics/lines/
ln -s /orange/adamginsburg/ACES/mosaics/*max_mosaic*png /orange/adamginsburg/web/secure/ACES/mosaics/lines/
"""
import shutil
import glob
import radio_beam
import os
from astropy import units as u
from astropy.io import fits
from astropy import log
from astropy.wcs import WCS
from aces.imaging.make_mosaic import make_mosaic, read_as_2d, get_peak, get_m0, all_lines
# import os
# from functools import partial
from multiprocessing import Process, Pool
import multiprocessing
import numpy as np
from spectral_cube.cube_utils import mosaic_cubes
import dask
from functools import partial
import multiprocessing.pool as mpp
from itertools import repeat

from aces import conf
import uvcombine

# from dask.distributed import Client
# client = Client(processes=True, n_workers=1, threads_per_worker=16)
# print(client)
# dask.config.set(scheduler=client)
if os.getenv('NO_PROGRESSBAR') is None and not (os.getenv('ENVIRON') == 'BATCH'):
    from dask.diagnostics import ProgressBar
    pbar = ProgressBar(minimum=20)  # don't show pbar <5s
    pbar.register()

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


def main():

    np.seterr('ignore')

    header = fits.Header.fromtextfile(f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr')

    processes = []
    for func in (residuals, reimaged, reimaged_high, continuum, cs21, hcop, hnco,
                 h40a):
        print(f"Starting function {func}")
        proc = Process(target=func, args=(header,))
        proc.start()
        processes.append(proc)

    for proc in processes:
        proc.join()

    # do this _after_
    all_lines(header)


def main_():

    np.seterr('ignore')

    header = fits.Header.fromtextfile(f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr')

    residuals(header)
    reimaged(header)
    reimaged_high(header)
    continuum(header)
    hcop(header)
    hnco(header)
    h40a(header)
    all_lines(header)


def continuum(header):
    log.info("12m continuum")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*cont.I.tt0.pbcor.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*25_27_29_31_33_35*cont.I.tt0.pbcor.fits')
    hdus = [read_as_2d(fn) for fn in filelist]
    print(flush=True)
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*I.pb.tt0.fits')
    weightfiles += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual.pb.tt0.fits')
    assert len(weightfiles) == len(filelist)
    wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='continuum_commonbeam_circular',
                commonbeam='circular',
                weights=wthdus,
                cbar_unit='Jy/beam', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.01, min_cut=-0.001),
                target_header=header,
                )
    print(flush=True)
    make_mosaic(hdus, name='continuum', weights=wthdus,
                cbar_unit='Jy/beam', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.01, min_cut=-0.001),
                target_header=header,
                )

    feath = uvcombine.feather_simple(f'{basepath}/mosaics/12m_continuum_commonbeam_circular_mosaic.fits',
                                     f'{basepath}/mosaics/7m_continuum_commonbeam_circular_mosaic.fits')
    fits.PrimaryHDU(data=feath.real,
                    header=fits.get_header(f'{basepath}/mosaics/12m_continuum_commonbeam_circular_mosaic.fits')
                    ).writeto(f'{basepath}/mosaics/feather_7m12m_continuum_commonbeam_circular_mosaic.fits',
                              overwrite=True)


def reimaged(header):
    log.info("12m continuum reimaged")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*25_27_29_31_33_35*cont.I.iter1.image.tt0.pbcor')
    hdus = [read_as_2d(fn) for fn in filelist]
    print(flush=True)
    weightfiles = [x.replace(".image.tt0.pbcor", ".pb.tt0") for x in filelist]
    weightfiles_ = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*25_27_29_31_33_35*I.iter1.pb.tt0')
    assert len(weightfiles) == len(filelist)
    for missing in set(weightfiles_) - set(weightfiles):
        print(f"Missing {missing}")
    wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='continuum_commonbeam_circular_reimaged',
                commonbeam='circular', weights=wthdus, cbar_unit='Jy/beam',
                array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.01,
                                 min_cut=-0.001),
                target_header=header,
                )
    print(flush=True)
    make_mosaic(hdus, name='continuum_reimaged', weights=wthdus,
                cbar_unit='Jy/beam', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.01, min_cut=-0.001),
                target_header=header,
                )


def reimaged_high(header):
    log.info("12m continuum reimaged")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33_35*cont.I.iter1.image.tt0.pbcor')
    hdus = [read_as_2d(fn) for fn in filelist]
    print(flush=True)
    weightfiles = [x.replace(".image.tt0.pbcor", ".pb.tt0") for x in filelist]
    weightfiles_ = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33_35*I.iter1.pb.tt0')
    assert len(weightfiles) == len(filelist)
    for missing in set(weightfiles_) - set(weightfiles):
        print(f"Missing {missing}")
    wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='continuum_commonbeam_circular_reimaged_spw33_35',
                commonbeam='circular', weights=wthdus, cbar_unit='Jy/beam',
                array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.01,
                                 min_cut=-0.001),
                target_header=header,
                )
    print(flush=True)
    make_mosaic(hdus, name='continuum_reimaged_spw33_35', weights=wthdus,
                cbar_unit='Jy/beam', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.01, min_cut=-0.001),
                target_header=header,
                )


def residuals(header):
    log.info("12m continuum residuals")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*25_27_29_31_33_35*cont.I.iter1.residual.tt0')
    hdus = [read_as_2d(fn) for fn in filelist]
    print(flush=True)
    weightfiles = [x.replace(".residual.tt0", ".pb.tt0") for x in filelist]
    beamfiles = [x.replace(".residual.tt0", ".image.tt0") for x in filelist]
    beams = radio_beam.Beams(beams=[radio_beam.Beam.from_fits_header(read_as_2d(fn)[0].header)
                                    for fn in beamfiles])
    assert len(weightfiles) == len(filelist)
    wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='continuum_residual_reimaged', weights=wthdus,
                cbar_unit='Jy/beam', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='linear', max_cut=0.001, min_cut=-0.001),
                target_header=header,
                )
    print(flush=True)
    cb = radio_beam.Beam.from_fits_header(fits.getheader(f'{basepath}/mosaics/12m_continuum_commonbeam_circular_mosaic.fits'))
    make_mosaic(hdus, name='continuum_residual_commonbeam_circular_reimaged',
                commonbeam=cb,
                beams=beams,
                weights=wthdus,
                cbar_unit='Jy/beam', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.001, min_cut=-0.001),
                target_header=header,
                )


def hcop(header):
    log.info("12m HCO+")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw29.cube.I.pbcor.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw29.cube.I.pbcor.fits')
    hdus = [get_peak(fn).hdu for fn in filelist]
    print(flush=True)
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*.Sgr_A_star_sci.spw29.mfs.I.pb.fits.gz')
    weightfiles += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*.Sgr_A_star_sci.spw29.mfs.I.pb.fits.gz')
    wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='hcop_max', cbar_unit='K', array='12m', basepath=basepath,
                weights=wthdus,
                norm_kwargs=dict(max_cut=20, min_cut=-0.5, ),
                target_header=header,
                )
    hdus = [get_m0(fn, suffix='_hcop').hdu for fn in filelist]
    print(flush=True)
    make_mosaic(hdus, name='hcop_m0', cbar_unit='K km/s', array='12m', basepath=basepath,
                weights=wthdus,
                target_header=header,
                norm_kwargs=dict(max_cut=100, min_cut=-10,))


def hnco(header):
    log.info("12m HNCO")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw31.cube.I.pbcor.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw31.cube.I.pbcor.fits')
    hdus = [get_peak(fn, suffix='_hnco').hdu for fn in filelist]
    print(flush=True)
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*.Sgr_A_star_sci.spw31.mfs.I.pb.fits.gz')
    weightfiles += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*.Sgr_A_star_sci.spw31.mfs.I.pb.fits.gz')
    wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='hnco_max', basepath=basepath, array='12m',
                target_header=header,
                weights=wthdus,
                norm_kwargs=dict(max_cut=20, min_cut=-0.5, ))
    hdus = [get_m0(fn, suffix='_hnco').hdu for fn in filelist]
    print(flush=True)
    make_mosaic(hdus, name='hnco_m0', cbar_unit='K km/s', array='12m',
                basepath=basepath, weights=wthdus, target_header=header,
                norm_kwargs=dict(max_cut=100, min_cut=-10, ))


def h40a(header):
    log.info("12m H40a")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw33.cube.I.pbcor.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw33.cube.I.pbcor.fits')
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.image.pbcor')
    hdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=99.02295 * u.GHz, suffix='_h40a').hdu for fn in filelist]
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.pb')
    wthdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=99.02295 * u.GHz, suffix='_h40a').hdu for fn in weightfiles]
    make_mosaic(hdus, name='h40a_max', cbar_unit='K',
                norm_kwargs=dict(max_cut=5, min_cut=-0.01, stretch='asinh'),
                array='12m', basepath=basepath, weights=wthdus, target_header=header)
    hdus = [get_m0(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=99.02295 * u.GHz).hdu for fn in filelist]
    make_mosaic(hdus, name='h40a_m0', cbar_unit='K km/s',
                norm_kwargs={'max_cut': 20, 'min_cut': -1, 'stretch': 'asinh'},
                array='12m', basepath=basepath, weights=wthdus, target_header=header)


def cs21(header):
    log.info("12m cs21")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw33.cube.I.pbcor.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw33.cube.I.pbcor.fits')
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.image.pbcor')
    hdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=97.98095330 * u.GHz, suffix='_cs21').hdu for fn in filelist]
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.pb')
    wthdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=97.98095330 * u.GHz, suffix='_cs21').hdu for fn in weightfiles]
    make_mosaic(hdus, name='cs21_max', cbar_unit='K',
                norm_kwargs=dict(max_cut=5, min_cut=-0.01, stretch='asinh'),
                array='12m', basepath=basepath, weights=wthdus, target_header=header)
    hdus = [get_m0(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=97.98095330 * u.GHz).hdu for fn in filelist]
    make_mosaic(hdus, name='cs21_m0', cbar_unit='K km/s',
                norm_kwargs={'max_cut': 20, 'min_cut': -1, 'stretch': 'asinh'},
                array='12m', basepath=basepath, weights=wthdus, target_header=header)


def istarmap(self, func, iterable, chunksize=1):
    """starmap-version of imap
    """
    self._check_running()
    if chunksize < 1:
        raise ValueError(
            "Chunksize must be 1+, not {0:n}".format(
                chunksize))

    task_batches = mpp.Pool._get_tasks(func, iterable, chunksize)
    result = mpp.IMapIterator(self)
    self._taskqueue.put(
        (
            self._guarded_task_generation(result._job,
                                          mpp.starmapstar,
                                          task_batches),
            result._set_length
        ))
    return (item for chunk in result for item in chunk)


mpp.Pool.istarmap = istarmap


def apply_kwargs(fn, kwargs):
    return fn(**kwargs)


def starstarmap_with_kwargs(pool, fn, kwargs_iter):
    args_for_starmap = zip(repeat(fn), kwargs_iter)
    return pool.istarmap(apply_kwargs, args_for_starmap)


def cs21_cube_mosaicing(test=False, verbose=True, channels='all'):
    """
    This takes too long as a full cube, so we have to do it slice-by-slice
    """
    from spectral_cube import SpectralCube
    from tqdm import tqdm

    print("Listing files", flush=True)
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33.cube.I.iter1.image.pbcor')
    if test:
        filelist = filelist[:16]
    header = fits.Header.fromtextfile(f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr')
    header['NAXIS'] = 3
    header['CDELT3'] = 1.4844932
    header['CUNIT3'] = 'km/s'
    header['CRVAL3'] = 0
    header['NAXIS3'] = 3 if test else int(np.ceil(300 / header['CDELT3']))
    header['CRPIX3'] = header['NAXIS3'] // 2
    header['CTYPE3'] = 'VRAD'
    header['RESTFRQ'] = 97.98095330e9
    header['SPECSYS'] = 'LSRK'

    print("Converting spectral units", flush=True)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cubes = [SpectralCube.read(fn,
                                   format='casa_image',
                                   use_dask=True).with_spectral_unit(u.km / u.s,
                                                                     velocity_convention='radio',
                                                                     rest_value=97.98095330 * u.GHz)
                 for fn in filelist]
        weightcubes = [(SpectralCube.read(fn.replace(".image.pbcor", ".pb"), format='casa_image', use_dask=True)
                        .with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=97.98095330 * u.GHz)
                        ) for fn in filelist]

    # flag out wild outliers
    # there are 2 as of writing
    print("Filtering out cubes with sketchy beams", flush=True)
    ok = [cube for cube in cubes if cube.beam.major < 2.4 * u.arcsec]
    cubes = [cube for k, cube in zip(ok, cubes) if k]
    weightcubes = [cube for k, cube in zip(ok, weightcubes) if k]

    beams = radio_beam.Beams(beams=[cube.beam for cube in cubes])
    commonbeam = beams.common_beam()
    header.update(commonbeam.to_header_keywords())

    ww = WCS(header)
    wws = ww.spectral

    if channels == 'all':
        channels = range(header['NAXIS3'])
    if True:
        pbar = tqdm(channels, desc='Channels (mosaic)') if verbose else channels
        for chan in pbar:
            theader = header.copy()
            theader['NAXIS3'] = 1
            theader['CRPIX3'] = 1
            theader['CRVAL3'] = wws.pixel_to_world(chan).to(u.km / u.s).value

            if verbose:
                pbar.set_description(f'Channel {chan}')

            chanfn = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/CS21_CubeMosaic_channel{chan}.fits'
            if os.path.exists(f'/orange/adamginsburg/ACES/mosaics/CS21_Channels/{os.path.basename(chanfn)}'):
                print(f"Skipping completed channel {chan}", flush=True)
            else:
                print(f"Starting mosaic_cubes for channel {chan}", flush=True)
                mosaic_cubes(cubes,
                             target_header=theader,
                             commonbeam=commonbeam,
                             weights=weightcubes,
                             spectral_block_size=None,
                             output_file=chanfn,
                             method='channel',
                             verbose=verbose
                             )
                shutil.move(chanfn, '/orange/adamginsburg/ACES/mosaics/CS21_Channels/')

    """
    # Failed experiment; MyPool needs to get restored for this to be considered
    if False:
        print("Assembling headers")
        # pool.map version:
        theaders = [header.copy() for ii in channels]
        for chan, theader in enumerate(theaders):
            theader.update({'NAXIS3': 1, 'CRPIX3': 1, 'CRVAL': wws.pixel_to_world(chan).to(u.km/u.s).value})

        print("Spinning up the pool")
        pool = MyPool(processes=int(os.getenv("SLURM_NPROCS")) if os.getenv('SLURM_NPROCS') else None)
        print(f"Pool has {pool._processes} jobs")
        #pool.map(partial(mosaic_cubes_kwargs, commonbeam=commonbeam, spectral_block_size=None, method='channel', verbose=False),
        #         [{'cubes':cubes, 'target_header': th, 'output_file': f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/CS21_CubeMosaic_channel{chan}.fits'}
        #          for chan, th in enumerate(theaders)]
        #         )
        rslt = starstarmap_with_kwargs(pool,
                                       partial(mosaic_cubes, commonbeam=commonbeam, spectral_block_size=None, method='channel', verbose=test),
                                       [{'cubes':cubes, 'target_header': th,
                                       'output_file': f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/CS21_CubeMosaic_channel{chan}.fits'}
                                       for chan, th in enumerate(theaders)]
                                      )
        for _ in tqdm(rslt, total=len(theaders), desc='starmap'):
            pass

        # Why use pool.imap then pool.wait?  for tqdm.
        pool.wait()
        pool.close()
        pool.join()
    """

    for chan in channels:
        chanfn = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/CS21_CubeMosaic_channel{chan}.fits'
        if os.path.exists(chanfn):
            shutil.move(chanfn, '/orange/adamginsburg/ACES/mosaics/CS21_Channels/')

    output_file = '/orange/adamginsburg/ACES/mosaics/CS21_CubeMosaic.fits'
    output_file = '/blue/adamginsburg/adamginsburg/ACES/workdir/CS21_CubeMosaic.fits'
    if not os.path.exists(output_file):
        # https://docs.astropy.org/en/stable/generated/examples/io/skip_create-large-fits.html#sphx-glr-generated-examples-io-skip-create-large-fits-py
        hdu = fits.PrimaryHDU(data=np.ones([5, 5, 5], dtype=float),
                              header=header,
                              )
        for kwd in ('NAXIS1', 'NAXIS2', 'NAXIS3'):
            hdu.header[kwd] = header[kwd]

        shape_opt = header['NAXIS3'], header['NAXIS2'], header['NAXIS1']

        header.tofile(output_file, overwrite=True)
        with open(output_file, 'rb+') as fobj:
            fobj.seek(len(header.tostring()) +
                      (np.prod(shape_opt) * np.abs(header['BITPIX'] // 8)) - 1)
            fobj.write(b'\0')

    hdu = fits.open(output_file, mode='update', overwrite=True)
    output_array = hdu[0].data
    hdu.flush()  # make sure the header gets written right

    pbar = tqdm(channels, desc='Channels')
    for chan in pbar:
        #chanfn = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/CS21_CubeMosaic_channel{chan}.fits'
        chanfn = f'/orange/adamginsburg/ACES/mosaics/CS21_Channels/CS21_CubeMosaic_channel{chan}.fits'
        if os.path.exists(chanfn):
            pbar.set_description('Channels (filling)')
            output_array[chan] = fits.getdata(chanfn)
            pbar.set_description('Channels (flushing)')
            hdu.flush()


def sio21_cube_mosaicing(test=False, verbose=True, channels='all'):
    """
    This takes too long as a full cube, so we have to do it slice-by-slice
    """
    from spectral_cube import SpectralCube
    from tqdm import tqdm

    print("Listing files", flush=True)
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw27.cube.I.iter1.image.pbcor')
    if test:
        filelist = filelist[:16]
    header = fits.Header.fromtextfile(f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr')
    header['NAXIS'] = 3
    header['CDELT3'] = 0.84
    header['CUNIT3'] = 'km/s'
    header['CRVAL3'] = 0
    header['NAXIS3'] = 3 if test else int(np.ceil(350 / header['CDELT3']))
    header['CRPIX3'] = header['NAXIS3'] // 2
    header['CTYPE3'] = 'VRAD'
    restfrq = 86.84696e9
    header['RESTFRQ'] = restfrq
    header['SPECSYS'] = 'LSRK'

    print("Converting spectral units", flush=True)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cubes = [SpectralCube.read(fn,
                                   format='casa_image',
                                   use_dask=True).with_spectral_unit(u.km / u.s,
                                                                     velocity_convention='radio',
                                                                     rest_value=restfrq * u.Hz)
                 for fn in filelist]
        weightcubes = [(SpectralCube.read(fn.replace(".image.pbcor", ".pb"), format='casa_image', use_dask=True)
                        .with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=restfrq * u.Hz)
                        ) for fn in filelist]

    # flag out wild outliers
    # there are 2 as of writing
    print("Filtering out cubes with sketchy beams", flush=True)
    ok = [cube for cube in cubes if cube.beam.major < 2.4 * u.arcsec]
    cubes = [cube for k, cube in zip(ok, cubes) if k]
    weightcubes = [cube for k, cube in zip(ok, weightcubes) if k]

    beams = radio_beam.Beams(beams=[cube.beam for cube in cubes])
    commonbeam = beams.common_beam()
    header.update(commonbeam.to_header_keywords())

    ww = WCS(header)
    wws = ww.spectral

    if channels == 'all':
        channels = range(header['NAXIS3'])
    if True:
        pbar = tqdm(channels, desc='Channels (mosaic)') if verbose else channels
        for chan in pbar:
            theader = header.copy()
            theader['NAXIS3'] = 1
            theader['CRPIX3'] = 1
            theader['CRVAL3'] = wws.pixel_to_world(chan).to(u.km / u.s).value

            if verbose:
                pbar.set_description(f'Channel {chan}')

            chanfn = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/SiO21_CubeMosaic_channel{chan}.fits'
            if os.path.exists(f'/orange/adamginsburg/ACES/mosaics/SiO21_Channels/{os.path.basename(chanfn)}'):
                print(f"Skipping completed channel {chan}", flush=True)
            else:
                print(f"Starting mosaic_cubes for channel {chan}", flush=True)
                mosaic_cubes(cubes,
                             target_header=theader,
                             commonbeam=commonbeam,
                             weights=weightcubes,
                             spectral_block_size=None,
                             output_file=chanfn,
                             method='channel',
                             verbose=verbose
                             )
                shutil.move(chanfn, '/orange/adamginsburg/ACES/mosaics/SiO21_Channels/')

    for chan in channels:
        chanfn = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/SiO21_CubeMosaic_channel{chan}.fits'
        if os.path.exists(chanfn):
            shutil.move(chanfn, '/orange/adamginsburg/ACES/mosaics/SiO21_Channels/')

    output_file = '/orange/adamginsburg/ACES/mosaics/SiO21_CubeMosaic.fits'
    output_file = '/blue/adamginsburg/adamginsburg/ACES/workdir/SiO21_CubeMosaic.fits'
    if not os.path.exists(output_file):
        # https://docs.astropy.org/en/stable/generated/examples/io/skip_create-large-fits.html#sphx-glr-generated-examples-io-skip-create-large-fits-py
        hdu = fits.PrimaryHDU(data=np.ones([5, 5, 5], dtype=float),
                              header=header,
                              )
        for kwd in ('NAXIS1', 'NAXIS2', 'NAXIS3'):
            hdu.header[kwd] = header[kwd]

        shape_opt = header['NAXIS3'], header['NAXIS2'], header['NAXIS1']

        header.tofile(output_file, overwrite=True)
        with open(output_file, 'rb+') as fobj:
            fobj.seek(len(header.tostring()) +
                      (np.prod(shape_opt) * np.abs(header['BITPIX'] // 8)) - 1)
            fobj.write(b'\0')

    hdu = fits.open(output_file, mode='update', overwrite=True)
    output_array = hdu[0].data
    hdu.flush()  # make sure the header gets written right

    pbar = tqdm(channels, desc='Channels')
    for chan in pbar:
        #chanfn = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/SiO21_CubeMosaic_channel{chan}.fits'
        chanfn = f'/orange/adamginsburg/ACES/mosaics/SiO21_Channels/SiO21_CubeMosaic_channel{chan}.fits'
        if os.path.exists(chanfn):
            pbar.set_description('Channels (filling)')
            output_array[chan] = fits.getdata(chanfn)
            pbar.set_description('Channels (flushing)')
            hdu.flush()


def hnco_cube_mosaicing(test=False, verbose=True, channels='all'):
    """
    This takes too long as a full cube, so we have to do it slice-by-slice
    """
    from spectral_cube import SpectralCube
    from tqdm import tqdm

    print("Listing files", flush=True)
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw31.cube.I.iter1.image.pbcor')
    if test:
        filelist = filelist[:16]
    header = fits.Header.fromtextfile(f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr')
    header['NAXIS'] = 3
    header['CDELT3'] = 0.21
    header['CUNIT3'] = 'km/s'
    header['CRVAL3'] = 0
    header['NAXIS3'] = 3 if test else int(np.ceil(250 / header['CDELT3']))
    header['CRPIX3'] = header['NAXIS3'] // 2
    header['CTYPE3'] = 'VRAD'
    restfrq = 87.925238e9
    header['RESTFRQ'] = restfrq
    header['SPECSYS'] = 'LSRK'

    print("Converting spectral units", flush=True)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cubes = [SpectralCube.read(fn,
                                   format='casa_image',
                                   use_dask=True).with_spectral_unit(u.km / u.s,
                                                                     velocity_convention='radio',
                                                                     rest_value=restfrq * u.Hz)
                 for fn in filelist]
        weightcubes = [(SpectralCube.read(fn.replace(".image.pbcor", ".pb"), format='casa_image', use_dask=True)
                        .with_spectral_unit(u.km / u.s, velocity_convention='radio', rest_value=restfrq * u.Hz)
                        ) for fn in filelist]

    # flag out wild outliers
    # there are 2 as of writing
    print("Filtering out cubes with sketchy beams", flush=True)
    ok = [cube for cube in cubes if cube.beam.major < 3.2 * u.arcsec]
    cubes = [cube for k, cube in zip(ok, cubes) if k]
    weightcubes = [cube for k, cube in zip(ok, weightcubes) if k]

    beams = radio_beam.Beams(beams=[cube.beam for cube in cubes])
    commonbeam = beams.common_beam()
    header.update(commonbeam.to_header_keywords())

    ww = WCS(header)
    wws = ww.spectral

    if channels == 'all':
        channels = range(header['NAXIS3'])
    if True:
        pbar = tqdm(channels, desc='Channels (mosaic)') if verbose else channels
        for chan in pbar:
            theader = header.copy()
            theader['NAXIS3'] = 1
            theader['CRPIX3'] = 1
            theader['CRVAL3'] = wws.pixel_to_world(chan).to(u.km / u.s).value

            if verbose:
                pbar.set_description(f'Channel {chan}')

            chanfn = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/HNCO_CubeMosaic_channel{chan}.fits'
            if os.path.exists(f'/orange/adamginsburg/ACES/mosaics/HNCO_Channels/{os.path.basename(chanfn)}'):
                print(f"Skipping completed channel {chan}", flush=True)
            else:
                print(f"Starting mosaic_cubes for channel {chan}", flush=True)
                mosaic_cubes(cubes,
                             target_header=theader,
                             commonbeam=commonbeam,
                             weights=weightcubes,
                             spectral_block_size=None,
                             output_file=chanfn,
                             method='channel',
                             verbose=verbose
                             )
                shutil.move(chanfn, '/orange/adamginsburg/ACES/mosaics/HNCO_Channels/')

    for chan in channels:
        chanfn = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/HNCO_CubeMosaic_channel{chan}.fits'
        if os.path.exists(chanfn):
            shutil.move(chanfn, '/orange/adamginsburg/ACES/mosaics/HNCO_Channels/')

    output_file = '/orange/adamginsburg/ACES/mosaics/HNCO_CubeMosaic.fits'
    output_file = '/blue/adamginsburg/adamginsburg/ACES/workdir/HNCO_CubeMosaic.fits'
    if not os.path.exists(output_file):
        # https://docs.astropy.org/en/stable/generated/examples/io/skip_create-large-fits.html#sphx-glr-generated-examples-io-skip-create-large-fits-py
        hdu = fits.PrimaryHDU(data=np.ones([5, 5, 5], dtype=float),
                              header=header,
                              )
        for kwd in ('NAXIS1', 'NAXIS2', 'NAXIS3'):
            hdu.header[kwd] = header[kwd]

        shape_opt = header['NAXIS3'], header['NAXIS2'], header['NAXIS1']

        header.tofile(output_file, overwrite=True)
        with open(output_file, 'rb+') as fobj:
            fobj.seek(len(header.tostring()) +
                      (np.prod(shape_opt) * np.abs(header['BITPIX'] // 8)) - 1)
            fobj.write(b'\0')

    hdu = fits.open(output_file, mode='update', overwrite=True)
    output_array = hdu[0].data
    hdu.flush()  # make sure the header gets written right

    pbar = tqdm(channels, desc='Channels')
    for chan in pbar:
        #chanfn = f'/blue/adamginsburg/adamginsburg/ACES/workdir/mosaics/HNCO_CubeMosaic_channel{chan}.fits'
        chanfn = f'/orange/adamginsburg/ACES/mosaics/HNCO_Channels/HNCO_CubeMosaic_channel{chan}.fits'
        if os.path.exists(chanfn):
            pbar.set_description('Channels (filling)')
            output_array[chan] = fits.getdata(chanfn)
            pbar.set_description('Channels (flushing)')
            hdu.flush()
