import numpy as np
import regions
from spectral_cube.spectral_cube import _regionlist_to_single_region
from astropy.table import Table
from astropy import units as u
from astropy.io import fits
from astropy import coordinates
from astropy import wcs
from astropy import log
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from make_mosaic import make_mosaic, read_as_2d, get_peak, get_m0

basepath = '/orange/adamginsburg/ACES/'

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

import glob

if __name__ == "__main__":

    log.info("12m continuum")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*cont.I.tt0.pbcor.fits')
    hdus = [read_as_2d(fn) for fn in filelist]
    print(flush=True)
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*25_27_29_31_33_35*I.pb.tt0.fits')
    assert len(weightfiles) == len(filelist)
    wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='continuum', weights=wthdus,
            cb_unit='Jy/beam', array='12m', basepath=basepath,
            norm_kwargs=dict(stretch='asinh', max_cut=0.01, min_cut=-0.001),
            )

    log.info("12m HCO+")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw29.cube.I.pbcor.fits')
    hdus = [get_peak(fn).hdu for fn in filelist]
    print(flush=True)
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*.Sgr_A_star_sci.spw29.mfs.I.pb.fits.gz')
    wthdus = [read_as_2d(fn, minval=0.3) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='hcop_max', cb_unit='K', array='12m', basepath=basepath,
                weights=wthdus,
                norm_kwargs=dict(max_cut=20, min_cut=-0.5, ),
               )
    hdus = [get_m0(fn).hdu for fn in filelist]
    print(flush=True)
    make_mosaic(hdus,  name='hcop_m0', cb_unit='K km/s', array='12m', basepath=basepath,
                weights=wthdus,
                norm_kwargs=dict(max_cut=100, min_cut=-10,))

    log.info("12m HNCO")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw31.cube.I.pbcor.fits')
    hdus = [get_peak(fn).hdu for fn in filelist]
    print(flush=True)
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*.Sgr_A_star_sci.spw31.mfs.I.pb.fits.gz')
    wthdus = [read_as_2d(fn, minval=0.3) for fn in weightfiles]
    print(flush=True)
    make_mosaic(hdus, name='hnco_max', basepath=basepath, array='12m',
                weights=wthdus,
                norm_kwargs=dict(max_cut=10, min_cut=-0.5, ))
    hdus = [get_m0(fn).hdu for fn in filelist]
    print(flush=True)
    make_mosaic(hdus, name='hnco_m0', cb_unit='K km/s', array='12m', basepath=basepath,
                weights=wthdus,
                norm_kwargs=dict(max_cut=100, min_cut=-10, ))

    log.info("12m H40a")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw33.cube.I.pbcor.fits')
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.image.pbcor')
    hdus = [get_peak(fn, slab_kwargs={'lo':-200*u.km/u.s, 'hi':200*u.km/u.s}, rest_value=99.02295*u.GHz).hdu for fn in filelist]
    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/calibrated/working//*spw33.cube.I.iter1.pb')
    wthdus = [get_peak(fn, slab_kwargs={'lo':-200*u.km/u.s, 'hi':200*u.km/u.s}, rest_value=99.02295*u.GHz).hdu for fn in weightfiles]
    make_mosaic(hdus, name='h40a_max', cb_unit='K', norm_kwargs=dict(max_cut=0.5, min_cut=-0.01, stretch='asinh'), array='12m', basepath=basepath, weights=wthdus)
    hdus = [get_m0(fn, slab_kwargs={'lo':-200*u.km/u.s, 'hi':200*u.km/u.s}, rest_value=99.02295*u.GHz).hdu for fn in filelist]
    make_mosaic(hdus, name='h40a_m0', cb_unit='K km/s', norm_kwargs={'max_cut': 20, 'min_cut':-1, 'stretch':'asinh'}, array='12m', basepath=basepath, weights=wthdus)

