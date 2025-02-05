import os
import glob
from astropy import units as u
from astropy.io import fits
from aces.imaging.make_mosaic import read_as_2d, get_peak, get_m0, rms
from aces.imaging.make_mosaic import (make_mosaic as make_mosaic_,
                                      all_lines as all_lines_,
                                      make_giant_mosaic_cube,
                                      make_downsampled_cube)

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


def make_giant_mosaic_cube_hcop_TP7m(**kwargs):
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
