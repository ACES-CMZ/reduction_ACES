import numpy as np
import warnings
import json
from astropy.table import Table,Column
from astropy import units as u
from astropy import log
from astropy import wcs
from astropy.io import fits
from astropy.stats import mad_std
from radio_beam import Beam
from radio_beam.beam import NoBeamException
from spectral_cube import SpectralCube
from spectral_cube.utils import NoBeamError, BeamWarning, StokesWarning
import scipy
import scipy.signal
from scipy import ndimage
import regions
import os
import glob
from functools import reduce
import operator
import re
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ..pipeline_scripts.merge_tclean_commands import commands
from ..retrieval_scripts.mous_map import get_mous_to_sb_mapping
from .. import conf
basepath = conf.basepath

warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning, append=True)
warnings.filterwarnings('ignore', category=BeamWarning, append=True)
warnings.filterwarnings('ignore', category=StokesWarning, append=True)
np.seterr(all='ignore')

def get_requested_sens():
    # use this file's path
    requested_fn = os.path.join(os.path.dirname(__file__), 'requested.txt')
    if not os.path.exists(requested_fn):
        requested_fn = f'{basepath}/reduction/analysis/requested.txt'
    from astropy.io import ascii
    tbl = ascii.read(requested_fn, data_start=2)
    return tbl

def get_psf_secondpeak(fn, show_image=False, min_radial_extent=1.5*u.arcsec,
                       max_radial_extent=5*u.arcsec, max_npix_peak=100,
                       specslice=slice(0,1)):
    """ REDUNDANT with get_psf_secondpeak_old, but this one is better

    Process:
        1. Find the first minimum of the PSF by taking the radial profile within 50 pixels
        2. Take the integral of the PSF within that range
        3. Calculate the residual of the PSF minus the CASA-fitted Gaussian beam
        4. Integrate that to get the fraction of flux outside the synthesized
        beam in the main lobe of the dirty beam
        5. Find the peak and the location of the peak residual

    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cube = SpectralCube.read(fn,
                                 format='casa_image' if not fn.endswith('.fits') else 'fits')
    psfim = cube[specslice][0]

    # this is a PSF: it _must_ have a peak of 1
    assert psfim.max() == 1.

    pixscale = wcs.utils.proj_plane_pixel_scales(cube.wcs.celestial)[0] * u.deg

    center = np.unravel_index(np.argmax(psfim), psfim.shape)
    cy, cx = center

    # can't extract a sub-region bigger than the cube
    for dim in psfim.shape:
        if max_npix_peak >= dim // 2:
            max_npix_peak = dim // 2 - 1

    npix = max_npix_peak

    # cannot cut out if the image is too small
    assert cx - npix >= 0
    assert cy - npix >= 0

    cutout = psfim[cy-npix:cy+npix+1, cx-npix:cx+npix+1]
    psfim_cutout = cutout

    try:
        beam = cube.beam
    except (AttributeError,NoBeamError):
        try:
            # assume we've appropriately sliced to get a single beam above
            beam = cube.beams[0]
        except (AttributeError,NoBeamError):
            log.error(f"File {fn} is not a valid PSF cube")
            return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)


    shape = cutout.shape
    sy, sx = shape

    fullbeam = beam.as_kernel(pixscale, x_size=npix*2+1, y_size=npix*2+1,)
    # will break things below fullbeam = beam.as_kernel(pixscale, x_size=sx, y_size=sy)

    Y, X = np.mgrid[0:sy, 0:sx]

    center = np.unravel_index(np.argmax(cutout), cutout.shape)
    cy, cx = center

    # elliptical version...
    dy = (Y - cy)
    dx = (X - cx)
    costh = np.cos(beam.pa)
    sinth = np.sin(beam.pa)
    rmajmin = beam.minor / beam.major

    rr = ((dx * costh + dy * sinth)**2 / rmajmin**2 +
          (dx * sinth - dy * costh)**2 / 1**2)**0.5

    rbin = (rr).astype(int)

    # assume the PSF first minimum is within 100 pixels of center
    radial_mean = ndimage.mean(cutout**2, labels=rbin, index=np.arange(max_npix_peak))

    # find the first negative peak (approximately); we include anything
    # within this radius as part of the main beam
    first_min_ind = scipy.signal.find_peaks(-radial_mean)[0][0]

    view = (slice(cy-first_min_ind.astype('int'), cy+first_min_ind.astype('int')+1),
            slice(cx-first_min_ind.astype('int'), cx+first_min_ind.astype('int')+1))
    data = cutout[view].value
    bm = fullbeam.array[view]
    # the data and beam must be concentric
    # and there must be only one peak location
    # (these checks are to avoid the even-kernel issue in which the center
    # of the beam can have its flux spread over four pixels)
    assert np.argmax(data) == np.argmax(bm)
    assert (bm.max() == bm).sum() == 1

    bmfit_residual = data-bm/bm.max()
    radial_mask = rr[view] < first_min_ind

    # calculate epsilon, the ratio of the PSF integral out to the first null to the integral of the PSF
    # the integral of the PSF should be very close to 1, but we want to peak-normalize to match the dirty beam
    synthbeam_integral = (fullbeam.array/fullbeam.array.max()).sum()
    log.debug(f"Synthetic beam integral = {synthbeam_integral}")
    dirtybeam_integral = (data / data.max() * radial_mask).sum()
    log.debug(f"Dirty beam integral = {dirtybeam_integral}")
    epsilon = synthbeam_integral / dirtybeam_integral
    log.debug(f"epsilon = {epsilon}")

    psf_integral_firstpeak = (data * radial_mask).sum()
    psf_residual_integral = (bmfit_residual * radial_mask).sum()
    residual_peak = bmfit_residual.max()
    residual_peak_loc = rr[view].flat[bmfit_residual.argmax()]

    peakloc_as = (residual_peak_loc * pixscale).to(u.arcsec)

    # pl.figure(3).clf()
    # bmradmean = ndimage.mean((fullbeam.array/fullbeam.array.max())**2, labels=rbin, index=np.arange(100))
    # pl.plot(radial_mean)
    # pl.plot(bmradmean)
    # pl.figure(1)

    # this finds the second peak
    # (useful for display)
    outside_first_peak_mask = (rr > first_min_ind) & (fullbeam.array < 1e-5)
    first_sidelobe_ind = scipy.signal.find_peaks(radial_mean *
                          (np.arange(len(radial_mean)) > first_min_ind))[0][0]
    max_sidelobe = cutout[outside_first_peak_mask].max()
    max_sidelobe_loc = cutout[outside_first_peak_mask].argmax()
    r_max_sidelobe = rr[outside_first_peak_mask][max_sidelobe_loc]

    if show_image:
        import pylab as pl

        # decide how big to make the plot
        if r_max_sidelobe * pixscale < min_radial_extent:
            radial_extent = (min_radial_extent / pixscale).decompose().value
        else:
            radial_extent = r_max_sidelobe
        if radial_extent * pixscale > max_radial_extent:
            radial_extent = (max_radial_extent / pixscale).decompose().value

        log.info(f"radial extent = {radial_extent},  "
                 f"r_max_sidelobe = {r_max_sidelobe}, "
                 "********" if r_max_sidelobe >  radial_extent else ""
                 f"first_sidelobe_ind={first_sidelobe_ind}, "
                 f"first_min_ind = {first_min_ind}")

        bm2 = beam.as_kernel(pixscale, x_size=radial_extent.astype('int')*2+1,
                             y_size=radial_extent.astype('int')*2+1,)
        view = (slice(cy-radial_extent.astype('int'), cy+radial_extent.astype('int')+1),
                slice(cx-radial_extent.astype('int'), cx+radial_extent.astype('int')+1))
        bmfit_residual2 = cutout[view].value-bm2.array/bm2.array.max()

        #extent = np.array([-first_min_ind, first_min_ind, -first_min_ind, first_min_ind])*pixscale.to(u.arcsec).value
        extent = np.array([-radial_extent, radial_extent, -radial_extent, radial_extent])*pixscale.to(u.arcsec).value
        ax = pl.gca()
        im = ax.imshow(bmfit_residual2, origin='lower',
                       interpolation='nearest', extent=extent, cmap='gray_r')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad="3%", axes_class=pl.matplotlib.axes.Axes)
        cb = pl.colorbar(mappable=im, cax=cax)
        pl.matplotlib.colorbar.ColorbarBase.add_lines(self=cb,
                                                      levels=[max_sidelobe],
                                                      colors=[(0.1,0.7,0.1,0.9)],
                                                      linewidths=1)

        ax.contour(bm2.array/bm2.array.max(), levels=[0.1,0.5,0.9], colors=['r']*3, extent=extent)
        ax.contour(rr[view], levels=[first_min_ind, r_max_sidelobe],
                   linestyles=['--',':'],
                   colors=[(0.2,0.2,1,0.5), (0.1,0.7,0.1,0.5)], extent=extent)
        ax.set_xlabel("RA Offset [arcsec]")
        ax.set_ylabel("Dec Offset [arcsec]")

    return (residual_peak,
            peakloc_as.value,
            psf_residual_integral/psf_integral_firstpeak,
            epsilon,
            first_min_ind*pixscale.to(u.arcsec),
            r_max_sidelobe*pixscale.to(u.arcsec),
            (rr, pixscale, cutout, beam, fullbeam, view, bmfit_residual)
           )






def imstats(fn, reg=None):
    """
    Compute a bunch of interesting image statistics and return them as a dictionary.

    Does more things if fn.replace('.image.tt0', '.psf.tt0') exists.

    Parameters
    ----------
    reg : regions.Region
        A region in which to calculate the noise (probably avoiding bright signal)

    """
    try:
        fh = fits.open(fn)
        data = fh[0].data
        ww = wcs.WCS(fh[0].header)
    except IsADirectoryError:
        cube = SpectralCube.read(fn, format='casa_image')
        data = cube[0].value
        ww = cube.wcs


    mad = mad_std(data, ignore_nan=True)
    peak = np.nanmax(data)
    imsum = np.nansum(data)
    sumgt5sig = np.nansum(data[data > 5*mad])
    sumgt3sig = np.nansum(data[data > 3*mad])

    pixscale = wcs.utils.proj_plane_pixel_area(ww)*u.deg**2

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)

        if 'cube' in locals():
            try:
                bm = cube.beam
                ppbeam = (bm.sr / pixscale).decompose()
                assert ppbeam.unit.is_equivalent(u.dimensionless_unscaled)
                ppbeam = ppbeam.value
            except NoBeamError:
                ppbeam = np.nan
                bm = Beam(np.nan)
        else:
            try:
                bm = Beam.from_fits_header(fh[0].header)
                ppbeam = (bm.sr / pixscale).decompose()
                assert ppbeam.unit.is_equivalent(u.dimensionless_unscaled)
                ppbeam = ppbeam.value
            except NoBeamException:
                ppbeam = np.nan
                bm = Beam(np.nan)


        meta = {'beam': bm.to_header_keywords(),
                'bmaj': bm.major.to(u.arcsec).value,
                'bmin': bm.minor.to(u.arcsec).value,
                'bpa': bm.pa.value,
                'mad': mad,
                'peak': peak,
                'peak/mad': peak / mad,
                'ppbeam': ppbeam,
                'sum': imsum,
                'fluxsum': imsum / ppbeam,
                'sumgt5sig': sumgt5sig,
                'sumgt3sig': sumgt3sig,
                'cellsize': (pixscale**0.5).to(u.arcsec).value
               }

    if reg is not None:
        try:
            reglist = regions.read_ds9(reg)
        except AttributeError: # newer version
            reglist = regions.Regions.read(reg)
        data = data.squeeze()
        composite_region = reduce(operator.or_, reglist)
        if hasattr(composite_region, 'to_mask'):
            msk = composite_region.to_mask()
        else:
            preg = composite_region.to_pixel(ww.celestial)
            msk = preg.to_mask()
        cutout_pixels = msk.cutout(data)[msk.data.astype('bool')]

        meta['mad_sample'] = mad_std(cutout_pixels, ignore_nan=True)
        meta['std_sample'] = np.nanstd(cutout_pixels)

    if fn.endswith('.image.tt0') or fn.endswith('.image.tt0.fits') or fn.endswith('.image.tt0.pbcor.fits') or fn.endswith('.image.tt0.pbcor'):
        psf_fn = fn.split(".image.tt0")[0] + ".psf.tt0"
    elif fn.endswith('.model.tt0') or fn.endswith('.model.tt0.fits') or fn.endswith('.model.tt0.pbcor.fits') or fn.endswith('.model.tt0.pbcor'):
        psf_fn = fn.split(".model.tt0")[0] + ".psf.tt0"
    elif fn.endswith('.image') or fn.endswith('.image.fits') or fn.endswith('.image.pbcor.fits') or fn.endswith('.image.pbcor'):
        psf_fn = fn.split(".image") + ".psf"
    else:
        raise IOError("Wrong image type passed to imstats: {fn}".format(fn=fn))

    if os.path.exists(psf_fn):
        try:
            psf_secondpeak, psf_secondpeak_loc, psf_sidelobe1_fraction, epsilon, firstmin, r_sidelobe, _ = get_psf_secondpeak(psf_fn)
        except IndexError:
            psf_secondpeak, psf_secondpeak_loc, psf_sidelobe1_fraction, epsilon, firstmin, r_sidelobe, _ = get_psf_secondpeak(psf_fn, max_npix_peak=200)
        meta['psf_secondpeak'] = psf_secondpeak
        meta['psf_epsilon'] = epsilon
        meta['psf_secondpeak_radius'] = psf_secondpeak_loc
        meta['psf_secondpeak_sidelobefraction'] = psf_sidelobe1_fraction
    else:
        meta['psf_secondpeak'] = np.nan
        meta['psf_epsilon'] = np.nan
        meta['psf_secondpeak_radius'] = np.nan
        meta['psf_secondpeak_sidelobefraction'] = np.nan

    return meta

"""
We want to calculate the "epsilon" value from https://ui.adsabs.harvard.edu/abs/1995AJ....110.2037J/abstract

epsilon = C1  / ( R2 - R1 )

C1 is the model convolved with the synthesized beam, summed
R1 is the residual, summed
R2 is the dirty map.  It can be calculated as R1 + model convolved with dirty beam.
epsion = C1 / D1
"""


def parse_fn(fn):
    """
    /orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15b4_X4f/calibrated/working/
    uid___A001_X15b4_X4f.s10_0.Sgr_A_star_sci.spw16_18_20_22_24_26.cont.I.iter1.image.tt0
    """    

    basename = os.path.basename(fn)

    split = basename.split(".")

    muid = split[0]
    mousmap = get_mous_to_sb_mapping('2021.1.00172.L')
    mousmap_ = {key.replace("/","_").replace(":","_"):val for key,val in mousmap.items()}
    sbname = mousmap_[muid]
    region = field = sbname.split("_")[3]

    robust = commands[sbname]['tclean_cont_pars']['aggregate']['robust']

    return {'region': region,
            'band': 'B3',
            'muid': muid,
            'array': '12Monly' if '12M' in split else '7M12M' if '7M12M' in split else '7M' if '7M' in split else '????',
            'robust': 'r'+str(robust),
            'suffix': split[-1],
            'bsens': 'bsens' in fn.lower(),
            'nobright': ('noco' in fn.lower()) or ('non2hp' in fn.lower()),
            'pbcor': 'pbcor' in fn.lower(),
           }

def assemble_stats(globstr, ditch_suffix=None):
    import glob
    from astropy.utils.console import ProgressBar

    allstats = []

    for fn in ProgressBar(glob.glob(globstr)):
        if fn.endswith('diff.fits') or fn.endswith('bsens-cleanest.fits'):
            continue
        if fn.count('.fits') > 1:
            # these are diff images, or something like that
            continue
        if ditch_suffix is not None:
            meta = parse_fn(fn.split(ditch_suffix)[0])
            # don't do this on the suffix-ditched version
            meta['pbcor'] = 'pbcor' in fn.lower()
        else:
            meta = parse_fn(fn)
        meta['filename'] = fn
        stats = imstats(fn, reg=get_noise_region(meta['region'], meta['band']))
        allstats.append({'meta': meta, 'stats': stats})

    return allstats


def get_noise_region(field, band):
    """
    Noise estimation regions should be in regions/noise_estimation_regions/{field}_noise_sampling.reg
    """
    try:
        basepath = os.path.dirname(__file__) + "/../"
        noisepath = os.path.join(basepath, 'noise_estimation_regions')
        assert os.path.exists(noisepath)
    except AssertionError:
        noisepath = f'{basepath}/reduction_ACES/aces/data/regions/noise_estimation_regions/'


    regfn = f"{noisepath}/{field}_noise_sampling.reg"

    if os.path.exists(regfn):
        return regfn



def get_psf_secondpeak_old(fn, neighborhood_size=5, threshold=0.01):

    from scipy import ndimage
    from scipy.ndimage import filters

    if fn.endswith('fits'):
        data = fits.getdata(fn)
    else:
        from casatools import image
        ia = image()
        ia.open(fn)
        data = ia.getchunk()
        ia.close()

    if data.ndim > 2:
        data = data.squeeze()
    if data.ndim > 2:
        data = data[0,:,:]

    data_max = filters.maximum_filter(data, neighborhood_size)

    maxima = (data == data_max)
    data_min = filters.minimum_filter(data, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    pkval = [data[slc].max() for slc in slices]

    if len(pkval) >= 2:
        secondmax = sorted(pkval)[-2]
        return secondmax
    else:
        return np.nan


class MyEncoder(json.JSONEncoder):
    "https://stackoverflow.com/a/27050186/814354"
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)


def savestats(basepath=basepath,
              suffix='image.tt0*', filetype=".fits"):
    
    stats = assemble_stats(f"{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/*/calibrated/working/*.cont.I.iter1.{suffix}{filetype}", ditch_suffix=f".{suffix[:-1]}")
    with open(f'{basepath}/tables/metadata_{suffix}.json', 'w') as fh:
        json.dump(stats, fh, cls=MyEncoder)

    requested = get_requested_sens()

    meta_keys = ['region', 'band', 'array', 'robust', 'suffix',
                 'bsens', 'pbcor', 'nobright', 'filename']
    stats_keys = ['bmaj', 'bmin', 'bpa', 'peak', 'sum', 'fluxsum', 'sumgt3sig',
                  'sumgt5sig', 'mad', 'mad_sample', 'std_sample', 'peak/mad',
                  'psf_secondpeak', 'psf_secondpeak_radius',
                  'psf_secondpeak_sidelobefraction', 'cellsize',
                 ]
    req_keys = ['B3_res', 'B3_sens', ]
    req_keys_head = ['Req_Res', 'Req_Sens']

    rows = []
    for entry in stats:
        band = entry['meta']['band']
        requested_this = requested[requested['Field'] == entry['meta']['region']]
        if len(requested_this) == 0:
            print(f"Skipped {entry['meta']['region']}")
            continue
        rows += [[entry['meta'][key] for key in meta_keys] +
                 [entry['stats'][key] if key in entry['stats'] else np.nan for key in stats_keys] +
                 [requested_this[key][0] for key in req_keys if band in key]
                ]

    tbl = Table(rows=rows, names=meta_keys+stats_keys+req_keys_head)

    # do some QA
    tbl.add_column(Column(name='SensVsReq', data=tbl['mad']*1e3/tbl['Req_Sens']))
    tbl.add_column(Column(name='BeamVsReq', data=(tbl['bmaj']*tbl['bmin'])**0.5/tbl['Req_Res']))
    tbl.add_column(Column(name='BmajVsReq', data=tbl['bmaj']/tbl['Req_Res']))

    tbl.write(f'{basepath}/tables/metadata_{suffix.strip("*")}.ecsv', overwrite=True)
    tbl.write(f'{basepath}/tables/metadata_{suffix.strip("*")}.html',
              format='ascii.html', overwrite=True)
    tbl.write(f'{basepath}/tables/metadata_{suffix.strip("*")}.tex', overwrite=True)
    tbl.write(f'{basepath}/tables/metadata_{suffix.strip("*")}.js.html',
              format='jsviewer')

    return tbl

def main():
    savestats()
