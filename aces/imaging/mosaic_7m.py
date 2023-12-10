import glob
from astropy import units as u
from astropy.io import fits
from aces.imaging.make_mosaic import make_mosaic, read_as_2d, get_peak, get_m0, all_lines, rms
from aces.imaging.make_mosaic import make_mosaic as make_mosaic_, all_lines as all_lines_

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

def make_mosaic(*args, folder='7m_flattened', **kwargs):
    return make_mosaic_(*args, folder=folder, **kwargs)


def all_lines(*args, folder='12m_flattened', **kwargs):
    return all_lines_(*args, folder=folder, **kwargs)


def main():

    header = fits.Header.fromtextfile(f'{basepath}/reduction_ACES/aces/imaging/data/header_7m.hdr')

    print("7m Continuum")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*16_18_20_22*cont.I.tt0.pbcor.fits')
    hdus = [read_as_2d(fn) for fn in filelist]

    weightfiles = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*16_18_20_22*I.pb.tt0.fits')
    assert len(weightfiles) == len(filelist)
    wthdus = [read_as_2d(fn, minval=0.5) for fn in weightfiles]
    print(flush=True)

    make_mosaic(hdus, name='continuum', norm_kwargs=dict(stretch='asinh',
                max_cut=0.2, min_cut=-0.025), cbar_unit='Jy/beam', array='7m',
                weights=wthdus,
                target_header=header,
                basepath=basepath,
                folder='continuum'
                )
    make_mosaic(hdus, name='continuum_commonbeam_circular',
                commonbeam='circular',
                weights=wthdus,
                cbar_unit='Jy/beam', array='7m', basepath=basepath,
                norm_kwargs=dict(stretch='asinh', max_cut=0.2, min_cut=-0.025),
                target_header=header,
                folder='continuum'
                )

    print("7m HCO+")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw20.cube.I.pbcor.fits')
    hdus = [get_peak(fn).hdu for fn in filelist]
    make_mosaic(hdus, name='hcop_max', cbar_unit='K', array='7m',
                target_header=header,
                basepath=basepath, norm_kwargs=dict(max_cut=5, min_cut=-0.1))
    hdus = [get_m0(fn).hdu for fn in filelist]
    make_mosaic(hdus, name='hcop_m0', cbar_unit='K km/s', array='7m',
                target_header=header,
                basepath=basepath, norm_kwargs=dict(min_cut=-25, max_cut=150))

    print("7m HNCO")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw22.cube.I.pbcor.fits')
    hdus = [get_peak(fn).hdu for fn in filelist]
    make_mosaic(hdus, name='hnco_max', array='7m', basepath=basepath,
                target_header=header,
                norm_kwargs=dict(max_cut=5, min_cut=-0.1))
    hdus = [get_m0(fn).hdu for fn in filelist]
    make_mosaic(hdus, name='hnco_m0', cbar_unit='K km/s', array='7m',
                target_header=header,
                basepath=basepath, norm_kwargs=dict(min_cut=-25, max_cut=150))

    print("7m H40a")
    filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*spw24.cube.I.pbcor.fits')
    hdus = [get_peak(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=99.02295 * u.GHz).hdu for fn in filelist]
    make_mosaic(hdus, name='h40a_max', cbar_unit='K',
                target_header=header,
                norm_kwargs=dict(max_cut=0.5, min_cut=-0.01, stretch='asinh'),
                array='7m', basepath=basepath)
    hdus = [get_m0(fn, slab_kwargs={'lo': -200 * u.km / u.s, 'hi': 200 * u.km / u.s}, rest_value=99.02295 * u.GHz).hdu for fn in filelist]
    make_mosaic(hdus, name='h40a_m0', cbar_unit='K km/s',
                target_header=header,
                norm_kwargs={'max_cut': 20, 'min_cut': -1, 'stretch': 'asinh'},
                array='7m', basepath=basepath)

    rms(prefix='7m_continuum', threshold=3)

    all_lines(header, array='7m', parallel=False)
