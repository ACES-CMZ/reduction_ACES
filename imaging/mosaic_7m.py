import reproject
from astropy.io import fits
from astropy import coordinates
from astropy import wcs
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd

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
filelist = glob.glob(f'{basepath}/rawdata/2021.1.00172.L/s*/g*/m*/product/*16_18_20_22*cont.I.tt0.pbcor.fits')

def read_as_2d(fn):
    fh = fits.open(fn)
    fh[0].data = fh[0].data.squeeze()
    ww = wcs.WCS(fh[0].header).celestial
    fh[0].header = ww.to_header()
    return fh

hdus = [read_as_2d(fn) for fn in filelist]

target_wcs, shape_out = find_optimal_celestial_wcs(hdus, frame='galactic')

array, footprint = reproject_and_coadd(hdus,
                                       target_wcs, shape_out=shape_out,
                                       reproject_function=reproject_interp)

fits.PrimaryHDU(data=array, header=target_wcs.to_header()).writeto(f'{basepath}/mosaics/7m_continuum_mosaic.fits', overwrite=True)

import pylab as pl
from astropy import visualization
pl.matplotlib.use('agg')
fig = pl.figure(figsize=(20,7))
ax = fig.add_subplot(111, projection=target_wcs)
im = ax.imshow(array, norm=visualization.simple_norm(array, stretch='asinh'))
pl.colorbar(mappable=im)
ax.coords[0].set_axislabel('Galactic Longitude')
ax.coords[1].set_axislabel('Galactic Latitude')

fig.savefig(f'{basepath}/mosaics/7m_continuum_mosaic.png', bbox_inches='tight')
