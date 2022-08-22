import glob
import radio_beam
from astropy import units as u
from astropy.io import fits
from astropy import log
from aces.imaging.make_mosaic import make_mosaic, read_as_2d, get_peak, get_m0
import numpy as np

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


def main():

    np.seterr('ignore')

    header = fits.Header.fromtextfile(f'{basepath}/reduction_ACES/aces/imaging/data/header_12m.hdr')

    residuals(header)
    reimaged(header)
    continuum(header)
    hcop(header)
    hnco(header)
    h40a(header)


def continuum(header):
    log.info("12m continuum")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*cont.I.tt0.pbcor.fits')
    hdus = [read_as_2d(fn) for fn in filelist]
    print(flush=True)
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*I.pb.tt0.fits')
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
                     commonbeam='circular',
                     weights=wthdus,
                     cbar_unit='Jy/beam', array='12m', basepath=basepath,
                     norm_kwargs=dict(stretch='asinh', max_cut=0.01, min_cut=-0.001),
                     target_header=header,
                     )
    print(flush=True)
    make_mosaic(hdus, name='continuum_reimaged', weights=wthdus,
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
                norm_kwargs=dict(stretch='asinh', max_cut=0.001, min_cut=-0.001),
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
    hdus = [get_peak(fn).hdu for fn in filelist]
    print(flush=True)
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*.Sgr_A_star_sci.spw29.mfs.I.pb.fits.gz')
    wthdus = [read_as_2d(fn, minval=0.3) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='hcop_max', cbar_unit='K', array='12m', basepath=basepath,
                weights=wthdus,
                norm_kwargs=dict(max_cut=20, min_cut=-0.5, ),
                target_header=header,
                )
    hdus = [get_m0(fn).hdu for fn in filelist]
    print(flush=True)
    make_mosaic(hdus, name='hcop_m0', cbar_unit='K km/s', array='12m', basepath=basepath,
                weights=wthdus,
                target_header=header,
                norm_kwargs=dict(max_cut=100, min_cut=-10,))


def hnco(header):
    log.info("12m HNCO")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw31.cube.I.pbcor.fits')
    hdus = [get_peak(fn).hdu for fn in filelist]
    print(flush=True)
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*.Sgr_A_star_sci.spw31.mfs.I.pb.fits.gz')
    wthdus = [read_as_2d(fn, minval=0.3) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='hnco_max', basepath=basepath, array='12m',
                target_header=header,
                weights=wthdus,
                norm_kwargs=dict(max_cut=10, min_cut=-0.5, ))
    hdus = [get_m0(fn).hdu for fn in filelist]
    print(flush=True)
    make_mosaic(hdus, name='hnco_m0', cbar_unit='K km/s', array='12m',
                basepath=basepath, weights=wthdus, target_header=header,
                norm_kwargs=dict(max_cut=100, min_cut=-10, ))


def h40a(header):
    log.info("12m H40a")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw33.cube.I.pbcor.fits')
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.image.pbcor')
    hdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=99.02295 * u.GHz).hdu for fn in filelist]
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.pb')
    wthdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=99.02295 * u.GHz).hdu for fn in weightfiles]
    make_mosaic(hdus, name='h40a_max', cbar_unit='K',
                norm_kwargs=dict(max_cut=0.5, min_cut=-0.01, stretch='asinh'),
                array='12m', basepath=basepath, weights=wthdus, target_header=header)
    hdus = [get_m0(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=99.02295 * u.GHz).hdu for fn in filelist]
    make_mosaic(hdus, name='h40a_m0', cbar_unit='K km/s',
                norm_kwargs={'max_cut': 20, 'min_cut': -1, 'stretch': 'asinh'},
                array='12m', basepath=basepath, weights=wthdus, target_header=header)
