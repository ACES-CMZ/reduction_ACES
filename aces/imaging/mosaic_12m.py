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
from aces.imaging.make_mosaic import (make_mosaic, read_as_2d, get_peak,
                                      get_m0, all_lines,
                                      make_downsampled_cube,
                                      make_giant_mosaic_cube,
                                      rms
                                      )
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


def logprint(x, **kwargs):
    print(x, flush=True, **kwargs)
    log.info(x)


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

    failure = False
    errors, tracebacks = [], []
    for proc in processes:
        proc.join()
        if proc.exitcode != 0:
            print(f"Exception caught from subprocess {proc}: exit code {proc.exitcode}")

    # do this _after_
    all_lines(header)

    if failure:
        print(errors)
        print(tracebacks)
        raise errors


def main_():

    np.seterr('ignore')

    header = fits.Header.fromtextfile(f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr')

    residuals(header)
    reimaged(header)
    reimaged_high(header)
    continuum(header)
    rms(prefix='12m_continuum', threshold=3)
    hcop(header)
    hnco(header)
    h40a(header)
    all_lines(header)


def continuum(header):
    logprint("12m continuum")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*cont.I.tt0.pbcor.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*25_27_29_31_33_35*.tt0.pbcor.fits')
    print("Read as 2d for files: ", end=None, flush=True)
    hdus = [read_as_2d(fn) for fn in filelist]
    for hdu, fn in zip(hdus, filelist):
        logprint(f'{fn}: {np.isnan(hdu.data).sum()}')
    logprint(filelist)
    print(flush=True)
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*I.pb.tt0.fits')
    weightfiles += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*.pb.tt0.fits')
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
                    header=fits.getheader(f'{basepath}/mosaics/12m_continuum_commonbeam_circular_mosaic.fits')
                    ).writeto(f'{basepath}/mosaics/feather_7m12m_continuum_commonbeam_circular_mosaic.fits',
                              overwrite=True)


def reimaged(header):
    logprint("12m continuum reimaged")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*25_27_29_31_33_35*cont.I*image.tt0.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*25_27_29_31_33_35*cont*tt0.pbcor.fits')
    print("Read as 2d for files: ", end=None, flush=True)
    hdus = [read_as_2d(fn) for fn in filelist]
    print(flush=True)
    logprint(filelist)
    weightfiles = [x.replace(".image.tt0.pbcor", ".pb.tt0") for x in filelist]
    weightfiles_ = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*25_27_29_31_33_35*I.iter1.pb.tt0')
    weightfiles_ += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*.pb.tt0.fits')
    assert len(weightfiles) == len(filelist)
    for missing in set(weightfiles_) - set(weightfiles):
        logprint(f"Missing {missing}")
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
    print(flush=True)
    make_mosaic(wthdus, name='primarybeam_coverage', weights=wthdus,
                cbar_unit='PB Level', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=1, min_cut=-0.001),
                target_header=header,
                )


def reimaged_high(header):
    logprint("12m continuum reimaged")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33_35*cont.I.iter1.image.tt0.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw33_35*cont*tt0.pbcor.fits')
    print("Read as 2d for files: ", end=None, flush=True)
    hdus = [read_as_2d(fn) for fn in filelist]
    print(flush=True)
    weightfiles = [x.replace(".image.tt0.pbcor", ".pb.tt0") for x in filelist]
    weightfiles_ = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33_35*I.iter1.pb.tt0')
    weightfiles_ += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw33_35*.pb.tt0.fits')
    assert len(weightfiles) == len(filelist)
    for missing in set(weightfiles_) - set(weightfiles):
        logprint(f"Missing {missing}")
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
    logprint("12m continuum residuals")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*25_27_29_31_33_35*cont.I.iter1.residual.tt0')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*25_27_29_31_33_35*cont.I*.residual.tt0')
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
    logprint("12m HCO+")
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
    logprint("12m HNCO")
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
    logprint("12m H40a")
    #filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw33.cube.I.pbcor.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw33.cube.I.pbcor.fits')
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw33.cube.I.pbcor.fits')
    hdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s},
                     rest_value=99.02295 * u.GHz,
                     suffix='_h40a').hdu for fn in filelist]
    for hdu, fn in zip(hdus, filelist):
        logprint(f'{fn}: {np.isnan(hdu.data).sum()}, {(hdu.data==0).sum()}')
    #weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.pb')
    weightfiles = [fn.replace(".iter1.image.pbcor", "iter1.pb") for fn in filelist]
    assert len(weightfiles) == len(filelist)
    wthdus = [get_peak(fn, slab_kwargs={'lo': -20 * u.km / u.s, 'hi': 20 * u.km / u.s},
                       rest_value=99.02295 * u.GHz,
                       threshold=0.5,
                       suffix='_h40a').hdu for fn in weightfiles]
    for hdu, fn in zip(wthdus, weightfiles):
        logprint(f'{fn}: {np.isnan(hdu.data).sum()}, {(hdu.data==0).sum()}')

    make_mosaic(hdus, name='h40a_max', cbar_unit='K',
                norm_kwargs=dict(max_cut=5, min_cut=-0.01, stretch='asinh'),
                array='12m', basepath=basepath, weights=wthdus, target_header=header)
    hdus = [get_m0(fn,
                   slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s},
                   rest_value=99.02295 * u.GHz).hdu for fn in filelist]
    for hdu, fn in zip(hdus, filelist):
        logprint(f'{fn}: {np.isnan(hdu.data).sum()}, {(hdu.data==0).sum()}')
    make_mosaic(hdus, name='h40a_m0', cbar_unit='K km/s',
                norm_kwargs={'max_cut': 20, 'min_cut': -1, 'stretch': 'asinh'},
                array='12m', basepath=basepath, weights=wthdus, target_header=header)


def cs21(header):
    logprint("12m cs21")
    #filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw33.cube.I.pbcor.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw33.cube.I.pbcor.fits')
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


def make_giant_mosaic_cube_cs21(**kwargs):
    """
    Sep 2023: Fields ar and ad are excluded because of their beams
    ad shouldn't be, but it is.
    """

    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33.cube.I.iter1.image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33.cube.I.manual*image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33.cube.I.iter1.reclean*image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*33.cube.I.manual.pbcor.fits')

    print(f"Found {len(filelist)} CS 2-1-containing spw33 files")

    restfrq = 97.98095330e9
    cdelt_kms = 1.4844932
    make_giant_mosaic_cube(filelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='CS21',
                           nchan=300,
                           beam_threshold=2.8 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/CS21_Channels/',
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        # for unknown reasons, this one _really_ doesn't like dask
        make_downsampled_cube(f'{basepath}/mosaics/cubes/CS21_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/CS21_CubeMosaic_downsampled9.fits',
                              overwrite=True,
                              use_dask=False
                              )


def make_giant_mosaic_cube_sio21(**kwargs):

    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw27.cube.I.iter1.image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw27.cube.I.manual*image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*27.cube.I.manual.pbcor.fits')

    print(f"Found {len(filelist)} SiO 2-1-containing spw27 files")

    restfrq = 86.84696e9
    cdelt_kms = 0.84455895
    make_giant_mosaic_cube(filelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='SiO21',
                           nchan=350,
                           beam_threshold=3.1 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/SiO21_Channels/',
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/SiO21_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/SiO21_CubeMosaic_downsampled9.fits',
                              )


def make_giant_mosaic_cube_hnco(**kwargs):

    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw31.cube.I.iter1.image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw31.cube.I.manual*image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*31.cube.I.manual.pbcor.fits')

    print(f"Found {len(filelist)} HNCO-containing spw31 files")

    restfrq = 87.925238e9
    cdelt_kms = 0.20818593  # smooth by 2 chans
    make_giant_mosaic_cube(filelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HNCO',
                           nchan=1200,
                           beam_threshold=3.2 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/HNCO_Channels/',
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HNCO_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HNCO_CubeMosaic_downsampled9.fits',
                              )


def make_giant_mosaic_cube_hc3n(**kwargs):

    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw35.cube.I.iter1.image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw35.cube.I.manual*image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*35.cube.I.manual.pbcor.fits')

    print(f"Found {len(filelist)} HC3N-containing spw35 files")

    restfrq = 100.0763e9
    cdelt_kms = 1.47015502
    make_giant_mosaic_cube(filelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HC3N',
                           nchan=300,
                           beam_threshold=2.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/HC3N_Channels/',
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HC3N_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HC3N_CubeMosaic_downsampled9.fits',
                              )


def make_giant_mosaic_cube_hnco_TP7m12m(**kwargs):

    filelist = sorted(glob.glob(f'{basepath}/upload/HNCO_comb_fits/12m_7m_TP_feather_cubes/Image_cubes/*.hnco43.image.fits'))
    weightfilelist = sorted(glob.glob(f'{basepath}/upload/HNCO_comb_fits/12m_7m_TP_feather_cubes/Weight_cubes/*.hnco43.image.weight.fits'))
    print(f"Found {len(filelist)} HNCO 7m+12m+TP FITS files")
    print(f"Found {len(weightfilelist)} HNCO 7m+12m+TP FITS weight files")
    assert len(weightfilelist) == len(filelist)
    for xx, yy in zip(filelist, weightfilelist):
        print(f'Beginning of filenames: {os.path.basename(xx.split(".")[0])}, {os.path.basename(yy.split(".")[0])}')
        assert os.path.basename(xx.split(".")[0]) == os.path.basename(yy.split(".")[0])

    restfrq = 87.925238e9
    #cdelt_kms = 0.10409296373
    cdelt_kms = 0.20818593  # smooth by 2 chans
    make_giant_mosaic_cube(filelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HNCO_7m12mTP',
                           nchan=1000,
                           beam_threshold=3.2 * u.arcsec,
                           target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m_bigpix.hdr',
                           channelmosaic_directory=f'{basepath}/mosaics/HNCO_7m12mTP_Channels/',
                           weightfilelist=weightfilelist,
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HNCO_7m12mTP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HNCO_7m12mTP_CubeMosaic_downsampled9.fits',
                              )


def make_giant_mosaic_cube_hcop(**kwargs):

    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw29.cube.I.iter1.image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw29.cube.I.manual*image.pbcor')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*29.cube.I.manual.pbcor.fits')

    print(f"Found {len(filelist)} HCOP-containing spw29 files")

    restfrq = 89.188526e9
    cdelt_kms = 0.20818593  # smooth by 2 chans
    make_giant_mosaic_cube(filelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HCOP',
                           nchan=1000,
                           beam_threshold=3.2 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/HCOP_Channels/',
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HCOP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HCOP_CubeMosaic_downsampled9.fits',
                              )


def make_giant_mosaic_cube_hnco_TP7m(**kwargs):

    raise NotImplementedError("Need weight cubes")
    filelist = sorted(glob.glob(f'{basepath}/upload/HNCO_comb_fits/7m_TP_feather_cubes/*.hnco43.image'))
    weightfilelist = sorted(glob.glob(f'{basepath}/upload/HNCO_comb_fits/7m_TP_feather_cubes/Weight_cubes/*.hnco43.image.weight.fits'))
    print(f"Found {len(filelist)} HNCO 7m+TP FITS files")
    print(f"Found {len(weightfilelist)} HNCO 7m+TP FITS weight files")
    assert len(weightfilelist) == len(filelist)
    for xx, yy in zip(filelist, weightfilelist):
        print(f'Beginning of filenames: {os.path.basename(xx.split(".")[0])}, {os.path.basename(yy.split(".")[0])}')
        assert os.path.basename(xx.split(".")[0]) == os.path.basename(yy.split(".")[0])

    restfrq = 87.925238e9
    #cdelt_kms = 0.10409296373
    cdelt_kms = 0.20818593  # smooth by 2 chans
    make_giant_mosaic_cube(filelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HNCO_7mTP',
                           nchan=1000,
                           beam_threshold=3.2 * u.arcsec,
                           target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_7m.hdr',
                           channelmosaic_directory=f'{basepath}/mosaics/HNCO_7mTP_Channels/',
                           weightfilelist=weightfilelist,
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HNCO_7mTP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HNCO_7mTP_CubeMosaic_downsampled9.fits',
                              )
