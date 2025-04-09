"""
After running this, symlink to web via:
ln -s /orange/adamginsburg/ACES/mosaics/12m_flattened/*m0_mosaic*png /orange/adamginsburg/web/secure/ACES/mosaics/lines/
ln -s /orange/adamginsburg/ACES/mosaics/12m_flattened/*max_mosaic*png /orange/adamginsburg/web/secure/ACES/mosaics/lines/
"""
import shutil
import glob
import radio_beam
import os
from astropy import units as u
from astropy.io import fits
from astropy import log
from astropy.wcs import WCS
from astropy.table import Table
from aces.imaging.make_mosaic import (read_as_2d, get_peak,
                                      get_m0,
                                      make_downsampled_cube,
                                      make_giant_mosaic_cube,
                                      rms,
                                      downsample_spectrally
                                      )
from aces.imaging.make_mosaic import make_mosaic as make_mosaic_, all_lines as all_lines_
# import os
# from functools import partial
from multiprocessing import Process, Pool
import multiprocessing
import numpy as np
import spectral_cube
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


def make_mosaic(*args, folder='12m_flattened', **kwargs):
    return make_mosaic_(*args, folder=folder, **kwargs)


def all_lines(*args, folder='12m_flattened', **kwargs):
    return all_lines_(*args, folder=folder, **kwargs)


def rms_(*args, **kwargs):
    return rms(prefix='12m_continuum',
               folder='continuum',
               threshold=3, **kwargs)


def main():

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--contonly", dest="contonly",
                      default=False,
                      action='store_true',
                      help="Run only continuum mosaicing?", metavar="contonly")
    parser.add_option("--remake-cache", dest="remake_cache",
                      default=False,
                      action='store_true',
                      help="Remake cached m0 and peak?", metavar="remake_cache")
    parser.add_option('--skip-cont', dest='skip_cont', default=False,
                      action='store_true',)
    (options, args) = parser.parse_args()

    np.seterr('ignore')

    header = fits.Header.fromtextfile(f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr')

    funcs = ((mosaic_field_am_pieces, residuals, reimaged, reimaged_high, continuum, rms_, manual_alpha)
             if options.contonly else
             (mosaic_field_am_pieces, residuals, reimaged, reimaged_high, continuum, rms_, manual_alpha, cs21, hcop, hnco, h40a))

    if not options.skip_cont:
        processes = []
        for func in funcs:
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
    if not options.contonly:
        print("Running all_lines")
        all_lines(header, use_cache=not options.remake_cache)

    if not options.skip_cont and failure:
        print(errors)
        print(tracebacks)
        raise errors


def main_():
    """
    I'm pretty sure main_ is the 'old version' that I was just commenting out without deleting...
    """

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


def check_files(filelist, funcname=None):
    if funcname is not None:
        logprint(f"Checking files for {funcname}")
    uidtb = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/aces_SB_uids.csv')
    for row in uidtb:
        matches = [(row['12m MOUS ID'] in fn) or
                   (f'_{row["Obs ID"]}.' in os.path.basename(fn))
                   for fn in filelist]
        print(row['Obs ID'], sum(matches))
        if sum(matches) != 1:
            for fn in filelist:
                if row['12m MOUS ID'] in fn:
                    print(fn)
            raise ValueError(f"Missing {row['Obs ID']} or too many (sum(matches)= {sum(matches)}), matches={[fn for fn in filelist if row['12m MOUS ID'] in fn]}")

    for fn in filelist:
        assert os.path.exists(fn
            .replace("image.pbcor.statcont.contsub.fits", "pb")
            .replace("image.pbcor", "pb")
        ), f"No pb found for {fn}"
        assert os.path.exists(fn
            .replace("image.pbcor.statcont.contsub.fits", "weight")
            .replace("image.pbcor", "weight")
        ), f"No weight found for {fn}"


def continuum(header):
    logprint("12m continuum: default/product version")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*cont.I.tt0.pbcor.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*cont.I.manual.pbcor.tt0.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*25_27_29_31_33_35*.tt0.pbcor.fits')

    check_files(filelist, funcname='continuum')

    print("CONTINUUM (default/product version) Reading as 2d for files: ", end=None, flush=True)
    hdus = [read_as_2d(fn) for fn in filelist]
    for hdu, fn in zip(hdus, filelist):
        if isinstance(hdu, fits.HDUList):
            hdu = hdu[0]
        logprint(f'{fn}: {np.isnan(hdu.data).sum()}')
    logprint(filelist)
    print(flush=True)
    #weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*I.pb.tt0.fits')
    #weightfiles += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*.pb.tt0.fits')
    weightfiles = [fn.replace(".image.tt0.pbcor", ".weight.tt0").replace(".I.tt0.pbcor", ".I.weight.tt0").replace('manual.pbcor.tt0', 'manual.weight.tt0')
                   for fn in filelist]
    # for product version, we need to use what we're given...
    weightfiles = [fn.replace(".weight.tt0.fits", ".pb.tt0.fits").replace(".weight.tt0", ".pb.tt0.fits") for fn in weightfiles if not os.path.exists(fn)]
    assert len(weightfiles) == len(filelist)
    wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='continuum_commonbeam_circular',
                commonbeam='circular',
                weights=wthdus,
                cbar_unit='Jy/beam', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.008, min_cut=-0.0002),
                target_header=header,
                folder='continuum'
                )
    print(flush=True)
    make_mosaic(hdus, name='continuum', weights=wthdus,
                cbar_unit='Jy/beam', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.008, min_cut=-0.0002),
                target_header=header,
                folder='continuum'
                )

    feath = uvcombine.feather_simple(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_mosaic.fits',
                                     f'{basepath}/mosaics/continuum/7m_continuum_commonbeam_circular_mosaic.fits')
    fits.PrimaryHDU(data=feath.real,
                    header=fits.getheader(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_mosaic.fits')
                    ).writeto(f'{basepath}/mosaics/continuum/feather_7m12m_continuum_commonbeam_circular_mosaic.fits',
                              overwrite=True)


def reimaged(header):
    logprint("12m continuum reimaged (glob is *.spw25_27_29_31_33_35.cont.I*image.tt0.pbcor)")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*.spw25_27_29_31_33_35.cont.I*image.tt0.pbcor')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*25_27_29_31_33_35*cont*tt0.pbcor.fits')
    
    # special-case field am
    am_index = [ii for ii, x in enumerate(filelist) if 'uid___A001_X15a0_X184' in x][0]
    #filelist[am_index] = f'{basepath}/mosaics/field_am/12m_field_am_mosaic.fits'
    filelist[am_index] = '/orange/adamginsburg/ACES/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X184/calibrated/working/uid___A001_X15a0_X184.Sgr_A_star_sci.cont.I.manual.lower_plus_upper.image.tt0.pbcor.fits'

    tt1filelist = [x.replace("tt0", "tt1") for x in filelist]

    check_files(filelist, funcname='reimaged')

    print("Reading as 2d for files (reimaged): ", end=None, flush=True)
    hdus = [read_as_2d(fn) for fn in filelist]
    print(flush=True)
    logprint(filelist)
    #weightfiles = [x.replace(".image.tt0.pbcor", ".pb.tt0") for x in filelist]
    #weightfiles_ = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*25_27_29_31_33_35*I.iter1.pb.tt0')
    #weightfiles_ += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*.pb.tt0.fits')
    weightfiles = [fn.replace(".image.tt0.pbcor", ".weight.tt0").replace(".I.tt0.pbcor", ".I.weight.tt0") for fn in filelist]
    pbfiles = [fn.replace(".image.tt0.pbcor", ".pb.tt0").replace(".I.tt0.pbcor", ".I.pb.tt0") for fn in filelist]
    #assert len(weightfiles) == len(filelist)
    #for missing in set(weightfiles_) - set(weightfiles):
    #    logprint(f"Missing {missing}")
    print("Reading as 2d for weightfiles (reimaged): ", end=None, flush=True)
    wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='continuum_commonbeam_circular_reimaged',
                commonbeam='circular', weights=wthdus, cbar_unit='Jy/beam',
                array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.008,
                                 min_cut=-0.0002),
                target_header=header,
                folder='continuum'
                )
    print(flush=True)
    make_mosaic(hdus, name='continuum_reimaged', weights=wthdus,
                cbar_unit='Jy/beam', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.008, min_cut=-0.0002),
                target_header=header,
                folder='continuum'
                )
    print(flush=True)
    make_mosaic(wthdus, name='average_weight',
                cbar_unit='Average Weight', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=1, min_cut=-0.0002),
                target_header=header,
                folder='continuum'
                )
    print(flush=True)
    print("Reading as 2d for pb (reimaged): ", end=None, flush=True)
    pbhdus = [read_as_2d(fn, minval=0.5) for fn in pbfiles]
    make_mosaic(pbhdus, name='primarybeam_coverage', weights=wthdus,
                cbar_unit='PB Level', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=1, min_cut=-0.0002),
                target_header=header,
                folder='continuum'
                )

    # feather with non-reimaged 7m (never did 7m reimaging)
    feath = uvcombine.feather_simple(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits',
                                     f'{basepath}/mosaics/continuum/7m_continuum_commonbeam_circular_mosaic.fits')
    fits.PrimaryHDU(data=feath.real,
                    header=fits.getheader(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits')
                    ).writeto(f'{basepath}/mosaics/continuum/feather_7m12m_continuum_commonbeam_circular_reimaged_mosaic.fits',
                              overwrite=True)

    # pbcorrect & mosaic tt1 image
    tt1hdus = [read_as_2d(fn) for fn in tt1filelist]
    pbhdus = [read_as_2d(fn) for fn in pbfiles]
    for tt1, pb in zip(tt1hdus, pbhdus):
        tt1[0].data /= pb[0].data
    make_mosaic(tt1hdus, name='continuum_commonbeam_circular_reimaged_tt1', weights=wthdus,
                cbar_unit='?', array='12m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=1, min_cut=-0.0002),
                target_header=header,
                folder='continuum',
                commonbeam='circular'
                )

    tt0 = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits')
    tt1 = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_tt1_mosaic.fits')
    alpha = tt1[0].data / tt0[0].data
    fits.PrimaryHDU(data=alpha, header=tt0[0].header).writeto(
        f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_alpha_mosaic.fits',
        overwrite=True)

    if os.path.exists(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_tt1_maskedrms_mosaic.fits'):
        ett0 = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_maskedrms_mosaic.fits')
        ett1 = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_tt1_maskedrms_mosaic.fits')
        ealpha = alpha * ((ett0[0].data / tt0[0].data)**2 + (ett1[0].data / tt1[0].data)**2)**0.5
        fits.PrimaryHDU(data=ealpha, header=tt0[0].header).writeto(
            f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_alphaerror_mosaic.fits',
            overwrite=True)


def reimaged_high(header, spws=('33_35', '25_27'), spw_names=('reimaged_high', 'reimaged_low')):
    for spw, name in zip(spws, spw_names):
        logprint(f"12m continuum {name} in reimaged_high")
        filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*.spw{spw}.cont.I*image.tt0.pbcor')
        filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*.spw{spw}.cont*tt0.pbcor.fits')

        check_files(filelist, funcname='reimaged_high')

        print(f"Reading as 2d for files (reimaged {name}): ", end=None, flush=True)
        hdus = [read_as_2d(fn) for fn in filelist]
        print(flush=True)
        #weightfiles = [x.replace(".image.tt0.pbcor", ".pb.tt0") for x in filelist]
        #weightfiles_ = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw{spw}*I.iter1.pb.tt0')
        #weightfiles_ += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw{spw}*.pb.tt0.fits')
        weightfiles = [fn.replace(".image.tt0.pbcor", ".weight.tt0").replace(".I.tt0.pbcor", ".I.weight.tt0") for fn in filelist]
        pbfiles = [fn.replace(".image.tt0.pbcor", ".pb.tt0").replace(".I.tt0.pbcor", ".I.pb.tt0") for fn in filelist]
        #assert len(weightfiles) == len(filelist)
        #for missing in set(weightfiles_) - set(weightfiles):
        #    logprint(f"Missing {missing}")
        wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
        print(flush=True)
        make_mosaic(hdus, name=f'continuum_commonbeam_circular_reimaged_spw{spw}',
                    commonbeam='circular', weights=wthdus, cbar_unit='Jy/beam',
                    array='12m', basepath=basepath,
                    norm_kwargs=dict(stretch='asinh', max_cut=0.008,
                                     min_cut=-0.0002),
                    target_header=header,
                    folder='continuum'
                    )
        print(flush=True)
        make_mosaic(hdus, name=f'continuum_reimaged_spw{spw}', weights=wthdus,
                    cbar_unit='Jy/beam', array='12m', basepath=basepath,
                    norm_kwargs=dict(stretch='asinh', max_cut=0.008, min_cut=-0.0002),
                    target_header=header,
                    folder='continuum'
                    )
        print(flush=True)
        make_mosaic(wthdus, name=f'average_weight_spw{spw}',
                    cbar_unit='Average Weight', array='12m', basepath=basepath,
                    norm_kwargs=dict(stretch='asinh', max_cut=1, min_cut=-0.0002),
                    target_header=header,
                    folder='continuum'
                    )
        print(flush=True)
        print("Reading as 2d for pb (reimaged): ", end=None, flush=True)
        pbhdus = [read_as_2d(fn, minval=0.5) for fn in pbfiles]
        make_mosaic(pbhdus, name=f'primarybeam_coverage_spw{spw}', weights=wthdus,
                    cbar_unit='PB Level', array='12m', basepath=basepath,
                    norm_kwargs=dict(stretch='asinh', max_cut=1, min_cut=-0.0002),
                    target_header=header,
                    folder='continuum'
                    )


def manual_alpha(header=None):
    """
    header is a dummy argument not used here, just added to allow manual_alpha to run
    """
    lofrqhdu = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_mosaic.fits')
    hifrqhdu = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits')
    lorms = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw25_27_maskedrms_mosaic.fits')
    lofrq = spectral_cube.Projection.from_hdu(lofrqhdu)
    hifrq = spectral_cube.Projection.from_hdu(hifrqhdu)

    hifrq_conv = hifrq.convolve_to(lofrq.beam)
    # convolving to lower resolution does not account for unit differences; both are in Jy/beam
    beam_area_ratio = lofrq.beam.sr / hifrq.beam.sr

    # frequencies quoted from paper table
    alpha = np.log10(hifrq_conv * beam_area_ratio / lofrq) / np.log10(99.53 / 86.54)

    fits.PrimaryHDU(data=alpha.value, header=lofrqhdu[0].header).writeto(
        f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_manual_alpha_mosaic.fits',
        overwrite=True)


def residuals(header):
    logprint("12m continuum residuals")
    for spw, name in zip(('25_27_29_31_33_35', '33_35', '25_27'), ('reimaged', 'reimaged_high', 'reimaged_low')):
        filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*.spw{spw}.*cont.I.iter1.residual.tt0')
        filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*.spw{spw}.*cont.I.manual.residual.tt0')
        filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*.spw{spw}.*cont.I*.residual.tt0')

        check_files(filelist, funcname='residuals')

        # check that field am, which is done, is included
        assert any([f'uid___A001_X15a0_X184.Sgr_A_star_sci.spw{spw}.cont.I.manual.residual.tt0' in fn
                    for fn in filelist])

        hdus = [read_as_2d(fn) for fn in filelist]
        print(flush=True)
        weightfiles = [x.replace(".residual.tt0", ".weight.tt0") for x in filelist]
        beamfiles = [x.replace(".residual.tt0", ".image.tt0") for x in filelist]
        beams = radio_beam.Beams(beams=[radio_beam.Beam.from_fits_header(read_as_2d(fn)[0].header)
                                        for fn in beamfiles])
        assert len(weightfiles) == len(filelist)
        wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
        print(flush=True)
        make_mosaic(hdus, name=f'continuum_residual_{name}', weights=wthdus,
                    cbar_unit='Jy/beam', array='12m', basepath=basepath,
                    norm_kwargs=dict(stretch='linear', max_cut=0.001, min_cut=-0.0002),
                    target_header=header,
                    )
        print(flush=True)
        if name == 'reimaged':
            cb = radio_beam.Beam.from_fits_header(fits.getheader(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits'))
        elif name == 'reimaged_high':
            cb = radio_beam.Beam.from_fits_header(fits.getheader(f'{basepath}/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_spw33_35_mosaic.fits'))
        else:
            raise NotImplementedError(f"name={name}")
        make_mosaic(hdus, name=f'continuum_residual_commonbeam_circular_{name}',
                    commonbeam=cb,
                    beams=beams,
                    weights=wthdus,
                    cbar_unit='Jy/beam', array='12m', basepath=basepath,
                    norm_kwargs=dict(stretch='asinh', max_cut=0.001, min_cut=-0.0002),
                    target_header=header,
                    )


def hcop(header):
    logprint("12m HCO+")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw29.cube.I.iter1*.image.pbcor.statcont.contsub.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw29*pbcor.statcont.contsub.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*sci29.cube.I*.image.pbcor.statcont.contsub.fits')

    check_files(filelist)

    hdus = [get_peak(fn).hdu for fn in filelist]
    print(flush=True)
    #weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*.Sgr_A_star_sci.spw29.mfs.I.pb.fits.gz')
    #weightfiles += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*.Sgr_A_star_sci.spw29.mfs.I.pb.fits.gz')
    weightfiles = [fn.replace("cube.I.pbcor.fits", "mfs.I.pb.fits.gz") for fn in filelist]
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
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw31.cube.I.iter1*.image.pbcor.statcont.contsub.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw31*pbcor.statcont.contsub.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*sci31.cube.I*.image.pbcor.statcont.contsub.fits')

    check_files(filelist)

    hdus = [get_peak(fn, suffix='_hnco').hdu for fn in filelist]
    print(flush=True)
    #weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*.Sgr_A_star_sci.spw31.mfs.I.pb.fits.gz')
    #weightfiles += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*.Sgr_A_star_sci.spw31.mfs.I.pb.fits.gz')
    weightfiles = [fn.replace("cube.I.pbcor.fits", "mfs.I.pb.fits.gz") for fn in filelist]
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
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33*.image.pbcor.statcont.contsub.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw33*pbcor.statcont.contsub.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*sci33.cube.I*.image.pbcor.statcont.contsub.fits')

    check_files(filelist)

    hdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s},
                     rest_value=99.02295 * u.GHz,
                     suffix='_h40a').hdu for fn in filelist]
    for hdu, fn in zip(hdus, filelist):
        logprint(f'{fn}: {np.isnan(hdu.data).sum()}, {(hdu.data == 0).sum()}')
    #weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.pb')
    weightfiles = [fn.replace(".iter1.image.pbcor", ".iter1.pb") for fn in filelist]
    assert len(weightfiles) == len(filelist)
    wthdus = [get_peak(fn, slab_kwargs={'lo': -20 * u.km / u.s, 'hi': 20 * u.km / u.s},
                       rest_value=99.02295 * u.GHz,
                       threshold=0.5,
                       suffix='_h40a').hdu for fn in weightfiles]
    for hdu, fn in zip(wthdus, weightfiles):
        logprint(f'{fn}: {np.isnan(hdu.data).sum()}, {(hdu.data == 0).sum()}')

    make_mosaic(hdus, name='h40a_max', cbar_unit='K',
                norm_kwargs=dict(max_cut=5, min_cut=-0.01, stretch='asinh'),
                array='12m', basepath=basepath, weights=wthdus, target_header=header)
    hdus = [get_m0(fn,
                   slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s},
                   rest_value=99.02295 * u.GHz).hdu for fn in filelist]
    for hdu, fn in zip(hdus, filelist):
        logprint(f'{fn}: {np.isnan(hdu.data).sum()}, {(hdu.data == 0).sum()}')
    make_mosaic(hdus, name='h40a_m0', cbar_unit='K km/s',
                norm_kwargs={'max_cut': 20, 'min_cut': -1, 'stretch': 'asinh'},
                array='12m', basepath=basepath, weights=wthdus, target_header=header)


def cs21(header):
    logprint("12m cs21")
    #filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw33.cube.I.pbcor.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw33.cube.I.pbcor.fits')
    #filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.image.pbcor')
    # can have 'reclean' files
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1*.image.pbcor.statcont.contsub.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*spw33*pbcor.statcont.contsub.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*sci33.cube.I*.image.pbcor.statcont.contsub.fits')

    check_files(filelist)

    hdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=97.98095330 * u.GHz, suffix='_cs21').hdu for fn in filelist]
    #weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.pb')
    weightfiles = [fn.replace("cube.I.iter1.image.pbcor", "cube.I.iter1.pb") for fn in filelist]
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


def get_weightfile(filename, spw):
    suffixes = ('.cube.I.iter1.reclean.weight.fits',
                '.cube.I.iter1.weight.fits',
                '.cube.I.manual.weight.fits',
                '.cube.I.manual.reclean.weight.fits',
                '.cube.I.iter1.reclean.weight',
                '.cube.I.iter1.weight',
                '.cube.I.manual.weight',
                '.cube.I.manual.reclean.weight',
               )
    spw_pre = ('.spw', '_sci')
    flist = []
    for pre in spw_pre:
        match_pre = False
        for suf in suffixes:
            matches = glob.glob(os.path.join(
                os.path.dirname(filename),
                f"*{pre}{spw}*{suf}"))
            if len(matches) > 0:
                flist += matches
                match_pre = True
                break
        if match_pre:
            break

    # dig through the sous/gous/mous for the appropriate weightfile
    if len(flist) != 1:
        uidtb = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/aces_SB_uids.csv')
        uidtb.add_index('Obs ID')
        sgmpath = os.path.join(basepath, 'data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_')
        if 'Sgr_A_st' in filename:
            fieldid = filename.split("Sgr_A_st_")[-1].split(".")[0]
            mousid = uidtb.loc[fieldid]['12m MOUS ID']
            for pre in spw_pre:
                match_pre = False
                for suf in suffixes:
                    matches = glob.glob(f'{sgmpath}{mousid}/calibrated/working/*{pre}{spw}*{suf}')
                    if len(matches) > 0:
                        flist += matches
                        match_pre = True
                        break
                if match_pre:
                    break

    # we have to have exactly 1 match
    assert len(flist) == 1, str(flist)
    return flist[0]


def make_giant_mosaic_cube_cs21(**kwargs):
    """
    Sep 2023: Fields ar and ad are excluded because of their beams
    ad shouldn't be, but it is.
    """

    #filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33.cube.I.iter1.image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*sci33.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33.cube.I.iter1.reclean*image.pbcor.statcont.contsub.fits')
    # this next line is only for field am and should be removed b/c we need the .pb/.weight
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*33.cube.I.manual.pbcor.fits')
    # TODO when feathered data are ready
    filelist = glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW33/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_33.image.statcont.contsub.fits')

    print(f"Found {len(filelist)} CS 2-1-containing spw33 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=33) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    # am has 2.83" beam
    restfrq = 97.98095330e9
    cdelt_kms = 1.4844932
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='CS21',
                           nchan=350,
                           beam_threshold=2.9 * u.arcsec,
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

    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw27.cube.I.iter1.image.pbcor.statcont.contsub.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw27.cube.I.manual*image.pbcor.statcont.contsub.fits')
    filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*sci27.cube.I.manual*image.pbcor.statcont.contsub.fits')
    # should exist in cal/working filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*27.cube.I.manual.pbcor.fits')
    filelist = glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW27/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_27.image.statcont.contsub.fits')

    print(f"Found {len(filelist)} SiO 2-1-containing spw27 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=27) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfrq = 86.84696e9
    cdelt_kms = 0.84455895
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='SiO21',
                           nchan=600,
                           # field am is 3.22"
                           beam_threshold=3.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/SiO21_Channels/',
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/SiO21_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/SiO21_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/SiO21_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/SiO21_CubeMosaic_spectrally.fits',
                              factor=2,
                              )


def make_giant_mosaic_cube_hnco(**kwargs):

    # filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw31.cube.I.iter1.image.pbcor.statcont.contsub.fits')
    # filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw31.cube.I.manual*image.pbcor.statcont.contsub.fits')
    # filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*sci31.cube.I.manual*image.pbcor.statcont.contsub.fits')
    # filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*31.cube.I.manual.pbcor.fits')
    filelist = glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW31/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.hnco43.image.statcont.contsub.fits')

    print(f"Found {len(filelist)} HNCO-containing spw31 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=31) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfrq = 87.925238e9
    cdelt_kms = 0.20818593  # smooth by 2 chans
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HNCO',
                           nchan=1400,
                           beam_threshold=3.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/HNCO_Channels/',
                           fail_if_cube_dropped=False,
                           #use_reproject_cube=True,
                           parallel=os.getenv('SLURM_NTASKS'),
                           **kwargs,)

    # 2024-06-15: using reproject_cube; no final step needed
    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HNCO_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HNCO_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/HNCO_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HNCO_CubeMosaic_spectrally.fits',
                              factor=7,
                              )


def make_giant_mosaic_cube_nsplus(**kwargs):

    filelist = glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW35/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_35.image.statcont.contsub.fits')

    print(f"Found {len(filelist)} NSplus-containing spw35 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=35) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfrq = 100.198e9
    cdelt_kms = 1.47015502
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='NSplus',
                           nchan=350,
                           # 2.35 is OK except o (2.38) and am (2.64)
                           # but am, despite its 2.64" beam, was excluded (??) with the 2.65" threshold.  As was o.  So maybe this is just a rerun?
                           beam_threshold=2.75 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/NSplus_Channels/',
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/NSplus_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/NSplus_CubeMosaic_downsampled9.fits',
                              )


def make_giant_mosaic_cube_hc3n(**kwargs):

    #filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw35.cube.I.iter1.image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw35.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*sci35.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw35.cube.I.iter1.reclean.image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*35.cube.I.manual.pbcor.fits')
    filelist = glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW35/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_35.image.statcont.contsub.fits')

    print(f"Found {len(filelist)} HC3N-containing spw35 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=35) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfrq = 100.0763e9
    cdelt_kms = 1.47015502
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HC3N',
                           nchan=350,
                           # 2.35 is OK except o (2.38) and am (2.64)
                           # but am, despite its 2.64" beam, was excluded (??) with the 2.65" threshold.  As was o.  So maybe this is just a rerun?
                           beam_threshold=2.75 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/HC3N_Channels/',
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HC3N_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HC3N_CubeMosaic_downsampled9.fits',
                              overwrite=True,
                              use_dask=False
                              )


def make_giant_mosaic_cube_hnco_TP7m12m_minitest(**kwargs):
    filelist = sorted(glob.glob(f'{basepath}/upload/Feather_12m_7m_TP/HNCO/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.hnco43.image.statcont.contsub.fits'))
    filelist = [x for x in filelist if any(y in x for y in ('_af.', '_s.', '_o.', '_f.'))]
    # what if we reverse order?
    filelist = filelist[::-1]
    print(f"Found {len(filelist)} HNCO 7m+12m+TP FITS files")
    log.debug(f"Full names of files: {filelist}")

    weightfilelist = [get_weightfile(fn, spw=31) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)
    print(f"Found {len(weightfilelist)} HNCO 7m+12m+TP FITS weight files")
    log.debug(f"Full names of weights: {weightfilelist}")
    assert len(weightfilelist) == len(filelist)
    for xx, yy in zip(filelist, weightfilelist):
        #print(f'Beginning of filenames: {os.path.basename(xx.split(".")[0])}, {os.path.basename(yy.split("/")[-1].split(".")[0])}')
        print(f"Basenames: {os.path.basename(xx)}, {os.path.basename(yy)}")

    restfrq = 87.925238e9
    #cdelt_kms = 0.10409296373
    cdelt_kms = 0.20818593  # smooth by 2 chans
    make_giant_mosaic_cube(filelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HNCO_7m12mTP_minitest',
                           nchan=400,
                           beam_threshold=3.3 * u.arcsec,
                           #target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m_bigpix_sofa.hdr',
                           channelmosaic_directory=f'{basepath}/mosaics/HNCO_7m12mTP_minitest_Channels/',
                           weightfilelist=weightfilelist,
                           #use_reproject_cube=True,
                           parallel=os.getenv('SLURM_NTASKS'),
                           fail_if_cube_dropped=False,
                           **kwargs,)


def make_giant_mosaic_cube_hnco_TP7m12m(**kwargs):

    #filelist = sorted(glob.glob(f'{basepath}/upload/HNCO_comb_fits/12m_7m_TP_feather_cubes/Image_cubes/*.hnco43.image.fits'))
    filelist = sorted(glob.glob(f'{basepath}/upload/Feather_12m_7m_TP/HNCO/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.hnco43.image.statcont.contsub.fits'))
    #weightfilelist = sorted(glob.glob(f'{basepath}/upload/HNCO_comb_fits/12m_7m_TP_feather_cubes/Weight_cubes/*.hnco43.image.weight.fits'))
    print(f"Found {len(filelist)} HNCO 7m+12m+TP FITS files")

    weightfilelist = [get_weightfile(fn, spw=31) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)
    print(f"Found {len(weightfilelist)} HNCO 7m+12m+TP FITS weight files")
    assert len(weightfilelist) == len(filelist)
    for xx, yy in zip(filelist, weightfilelist):
        print(f'Beginning of filenames: {os.path.basename(xx.split(".")[0])}, {os.path.basename(yy.split("/")[-1].split(".")[0])}')
        # These won't match b/c some are from Dan's upload, some are not
        #assert os.path.basename(xx.split(".")[0]) == os.path.basename(yy.split(".")[0])

    restfrq = 87.925238e9
    #cdelt_kms = 0.10409296373
    cdelt_kms = 0.20818593  # smooth by 2 chans
    make_giant_mosaic_cube(filelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HNCO_7m12mTP',
                           nchan=1400,
                           beam_threshold=3.3 * u.arcsec,
                           #target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_12m_bigpix.hdr',
                           channelmosaic_directory=f'{basepath}/mosaics/HNCO_7m12mTP_Channels/',
                           weightfilelist=weightfilelist,
                           #use_reproject_cube=True,
                           parallel=os.getenv('SLURM_NTASKS'),
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HNCO_7m12mTP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HNCO_7m12mTP_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/HNCO_7m12mTP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HNCO_7m12mTP_CubeMosaic_spectrally.fits',
                              factor=7,
                              )


def make_giant_mosaic_cube_hcop(**kwargs):

    #filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw29.cube.I.iter1.image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw29.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*sci29.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*29.cube.I.manual.pbcor.fits')
    filelist = glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW29/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.hco+10.image.statcont.contsub.fits')

    print(f"Found {len(filelist)} HCOP-containing spw29 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=29) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfrq = 89.188526e9
    cdelt_kms = 0.20818593  # smooth by 2 chans
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HCOP',
                           nchan=1400,
                           beam_threshold=3.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/HCOP_Channels/',
                           fail_if_cube_dropped=False,
                           #use_reproject_cube=True,
                           parallel=os.getenv('SLURM_NTASKS'),
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HCOP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HCOP_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/HCOP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HCOP_CubeMosaic_spectrally.fits',
                              factor=7,
                              )


def make_giant_mosaic_cube_hnco_TP7m(**kwargs):
    """
    Not 12m!
    """

    filelist = sorted(glob.glob(f'{basepath}/upload/HNCO_comb_fits/7m_TP_feather_cubes/*.hnco43.image'))
    filelist += sorted(glob.glob(f'{basepath}/upload/HNCO_comb_fits/7m_TP_feather_cubes/*.HNCO.image.fits'))
    filelist = sorted(filelist)

    weightfilelist = sorted(glob.glob(f'{basepath}/upload/HNCO_comb_fits/7m_TP_feather_cubes/HNCO_7M_weights/*.hnco43.image.weight.fits'))
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
                           nchan=1400,
                           beam_threshold=25 * u.arcsec,
                           target_header=f'{basepath}/reduction_ACES/aces/imaging/data/header_7m.hdr',
                           channelmosaic_directory=f'{basepath}/mosaics/HNCO_7mTP_Channels/',
                           weightfilelist=weightfilelist,
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HNCO_7mTP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HNCO_7mTP_CubeMosaic_downsampled9.fits',
                              overwrite=True
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/HNCO_7mTP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HNCO_7mTP_CubeMosaic_spectrally.fits',
                              factor=7,
                              )


def make_giant_mosaic_cube_ch3cho(**kwargs):

    filelist = sorted(glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW33/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_33.image.statcont.contsub.fits'))

    print(f"Found {len(filelist)} SO 3-2-containing spw33 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=33) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfreq = 98.900951e9
    cdelt_kms = 1.4844932
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfreq,
                           cdelt_kms=cdelt_kms,
                           cubename='CH3CHO',
                           nchan=350,
                           beam_threshold=3.1 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/CH3CHO_Channels/',
                           # we prefer to fail_if_dropped; enabling this for debugging
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/CH3CHO_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/CH3CHO_CubeMosaic_downsampled9.fits',
                              )


def make_giant_mosaic_cube_so32(**kwargs):

    #filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33.cube.I.iter1.image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*sci33.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw33.cube.I.iter1.reclean*image.pbcor.statcont.contsub.fits')
    ##filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*33.cube.I.manual.pbcor.fits')
    #filelist = sorted(filelist)
    filelist = sorted(glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW33/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_33.image.statcont.contsub.fits'))

    print(f"Found {len(filelist)} SO 3-2-containing spw33 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=33) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfreq = 99.29987e9
    cdelt_kms = 1.4844932
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfreq,
                           cdelt_kms=cdelt_kms,
                           cubename='SO32',
                           nchan=350,
                           beam_threshold=3.1 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/SO32_Channels/',
                           # we prefer to fail_if_dropped; enabling this for debugging
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/SO32_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/SO32_CubeMosaic_downsampled9.fits',
                              )


def make_giant_mosaic_cube_h13cn(**kwargs):

    filelist = sorted(glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW25/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_25.image.statcont.contsub.fits'))

    print(f"Found {len(filelist)} H13CN-containing spw25 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=25) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfreq = 86.33992e9
    cdelt_kms = 0.84455895
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfreq,
                           cdelt_kms=cdelt_kms,
                           cubename='H13CN',
                           nchan=600,
                           beam_threshold=3.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/H13CN_Channels/',
                           # we prefer to fail_if_dropped; enabling this for debugging
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/H13CN_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/H13CN_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/H13CN_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/H13CN_CubeMosaic_spectrally.fits',
                              factor=2,
                              )


def make_giant_mosaic_cube_h13cop(**kwargs):

    filelist = sorted(glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW27/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_27.image.statcont.contsub.fits'))

    print(f"Found {len(filelist)} h13cop-containing spw27 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=27) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfreq = 86.7543e9
    cdelt_kms = 0.84455895
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfreq,
                           cdelt_kms=cdelt_kms,
                           cubename='H13COp',
                           nchan=600,
                           beam_threshold=3.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/H13COp_Channels/',
                           # we prefer to fail_if_dropped; enabling this for debugging
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/H13COp_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/H13COp_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/H13COp_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/H13COp_CubeMosaic_spectrally.fits',
                              factor=2,
                              )


def make_giant_mosaic_cube_hn13c(**kwargs):

    filelist = sorted(glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW27/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_27.image.statcont.contsub.fits'))

    print(f"Found {len(filelist)} hn13c-containing spw27 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=27) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    # biggest beams are am & r (3.15")
    restfreq = 87.09085e9
    cdelt_kms = 0.84455895
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfreq,
                           cdelt_kms=cdelt_kms,
                           cubename='HN13C',
                           nchan=600,
                           beam_threshold=3.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/HN13C_Channels/',
                           # we prefer to fail_if_dropped; enabling this for debugging
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HN13C_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HN13C_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/HN13C_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HN13C_CubeMosaic_spectrally.fits',
                              factor=2,
                              )


def make_giant_mosaic_cube_so21(**kwargs):

    filelist = sorted(glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW25/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_25.image.statcont.contsub.fits'))

    print(f"Found {len(filelist)} SO 2-1-containing spw25 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=25) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfreq = 86.09395e9
    cdelt_kms = 0.84455895

    # biggest beams are am & r (3.15")
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfreq,
                           cdelt_kms=cdelt_kms,
                           cubename='SO21',
                           nchan=600,
                           beam_threshold=3.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/SO21_Channels/',
                           # we prefer to fail_if_dropped; enabling this for debugging
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/SO21_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/SO21_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/SO21_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/SO21_CubeMosaic_spectrally.fits',
                              factor=2,
                              )


def make_giant_mosaic_cube_hc15n(**kwargs):

    filelist = sorted(glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW25/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_25.image.statcont.contsub.fits'))

    print(f"Found {len(filelist)} HC15N-containing spw25 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=25) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfreq = 86.05496e9
    cdelt_kms = 0.84455895
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfreq,
                           cdelt_kms=cdelt_kms,
                           cubename='HC15N',
                           nchan=600,
                           beam_threshold=3.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/HC15N_Channels/',
                           # we prefer to fail_if_dropped; enabling this for debugging
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HC15N_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HC15N_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/HC15N_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HC15N_CubeMosaic_spectrally.fits',
                              factor=2,
                              )


def make_giant_mosaic_cube_h40a(**kwargs):

    filelist = sorted(glob.glob('/orange/adamginsburg/ACES/upload/Feather_12m_7m_TP/SPW33/cubes/Sgr_A_st_*.TP_7M_12M_feather_all.SPW_33.image.statcont.contsub.fits'))

    print(f"Found {len(filelist)} H40a-containing spw33 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=33) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfreq = 99.02295e9
    cdelt_kms = 1.4844932
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfreq,
                           cdelt_kms=cdelt_kms,
                           cubename='H40a',
                           nchan=600,
                           beam_threshold=3.1 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/H40a_Channels/',
                           # we prefer to fail_if_dropped; enabling this for debugging
                           fail_if_cube_dropped=False,
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/H40a_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/H40a_CubeMosaic_downsampled9.fits',
                              )


def make_giant_mosaic_cube_hcop_noTP(**kwargs):

    #filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw29.cube.I.iter1.image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw29.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*sci29.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*29.cube.I.manual.pbcor.fits')
    filelist = glob.glob('/orange/adamginsburg/ACES/upload/HCOp_feather_images/TM_7M_only_feather/Sgr_A_st_*.7M_12M_feather.SPW_29.image.statcont.contsub.fits')

    print(f"Found {len(filelist)} HCOP-containing spw29 files")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=29) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfrq = 89.188526e9
    cdelt_kms = 0.20818593  # smooth by 2 chans
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HCOP_noTP',
                           nchan=1400,
                           beam_threshold=3.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/HCOP_Channels/',
                           fail_if_cube_dropped=False,
                           #use_reproject_cube=True,
                           parallel=os.getenv('SLURM_NTASKS'),
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HCOP_noTP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HCOP_noTP_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/HCOP_noTP_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HCOP_noTP_CubeMosaic_spectrally.fits',
                              factor=7,
                              )


def make_giant_mosaic_cube_hcop_mopra(**kwargs):

    #filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw29.cube.I.iter1.image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*spw29.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working/*sci29.cube.I.manual*image.pbcor.statcont.contsub.fits')
    #filelist += glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/manual/*29.cube.I.manual.pbcor.fits')
    filelist = glob.glob('/orange/adamginsburg/ACES/upload/HCOp_feather_images/MOPRA_12M_7M_feather/Sgr_A_st_*.SD_7M_12M_feather.hco*.image.statcont.contsub.fits')

    print(f"Found {len(filelist)} HCOP-containing spw29 files (MOPRA)")

    check_files(filelist)

    weightfilelist = [get_weightfile(fn, spw=29) for fn in filelist]
    for fn in weightfilelist:
        assert os.path.exists(fn)

    restfrq = 89.188526e9
    cdelt_kms = 1.81607449832
    make_giant_mosaic_cube(filelist,
                           weightfilelist=weightfilelist,
                           reference_frequency=restfrq,
                           cdelt_kms=cdelt_kms,
                           cubename='HCOP_mopra',
                           nchan=165,
                           beam_threshold=3.3 * u.arcsec,
                           channelmosaic_directory=f'{basepath}/mosaics/HCOP_Channels/',
                           fail_if_cube_dropped=False,
                           #use_reproject_cube=True,
                           parallel=os.getenv('SLURM_NTASKS'),
                           **kwargs,)

    if not kwargs.get('skip_final_combination') and not kwargs.get('test'):
        make_downsampled_cube(f'{basepath}/mosaics/cubes/HCOP_mopra_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HCOP_mopra_CubeMosaic_downsampled9.fits',
                              )
        downsample_spectrally(f'{basepath}/mosaics/cubes/HCOP_mopra_CubeMosaic.fits',
                              f'{basepath}/mosaics/cubes/HCOP_mopra_CubeMosaic_spectrally.fits',
                              factor=3,
                              )


def mosaic_field_am_pieces():
    path = f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X184/calibrated/working/'

    basefn = 'Sgr_A_star_sci.spw25_27_29_31_33_35.cont.I.manual.{uplo}.{suffix}.fits'

    hdus = [read_as_2d(path + basefn.format(uplo=uplo, suffix='image.tt0.pbcor')) for uplo in ('lower', 'upper')]
    wthdus = [read_as_2d(path + basefn.format(uplo=uplo, suffix='weight.tt0')) for uplo in ('lower', 'upper')]
    mos = make_mosaic(hdus, name='field_am', array='12m', folder='field_am', basepath=conf.basepath, doplots=False, commonbeam=True)

    print(radio_beam.Beam.from_fits_header('/orange/adamginsburg/ACES//mosaics/field_am/12m_field_am_mosaic.fits'))
    os.link('/orange/adamginsburg/ACES//mosaics/field_am/12m_field_am_mosaic.fits',
            f'{path}/{basefn.format(uplo="lower_plus_upper", suffix="image.tt0.pbcor.fits")}')