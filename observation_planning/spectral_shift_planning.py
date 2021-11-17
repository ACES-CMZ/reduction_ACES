import numpy as np
from astropy import units as u
from spectral_cube import SpectralCube
from astropy.table import Table
import regions
from spectral_cube.spectral_cube import _regionlist_to_single_region


hcopcube = SpectralCube.read('/orange/adamginsburg/cmz/mopra/CMZ_3mm_HCO+.fits')
hncocube = SpectralCube.read('/orange/adamginsburg/cmz/mopra/CMZ_3mm_HNCO.fits')

tbl = Table.read('../SB_naming.tsv', format='ascii.csv', delimiter='\t')

vwidth = 102.585 * 2048 *u.m/u.s

mask = hcopcube.mask.include()
print(mask.sum())
bmask = np.zeros_like(mask)

composites = []

for row in tbl:
    indx = row['Proposal ID'][3:]
    print(f"Row {indx}: mask sum = {bmask.sum()}, which is {bmask.sum() / mask.sum() * 100:0.2f} percent of original")

    regs = regions.Regions.read(f'../regions/final_cmz{indx}.reg')
    pregs = [reg.to_pixel(hcopcube.wcs.celestial) for reg in regs]

    composite = _regionlist_to_single_region(pregs)
    composite.meta['label'] = indx
    composites.append(composite)

    voff = row['velocity offset']*u.km/u.s
    vmin, vmax = voff-vwidth/2, voff+vwidth/2

    insidespecslice = (hcopcube.spectral_axis > vmin) & (hcopcube.spectral_axis < vmax)
    outsidespecslice = (hcopcube.spectral_axis < vmin) & (hcopcube.spectral_axis > vmax)
    print(f"Specslice includes {insidespecslice.sum()} and excludes {outsidespecslice.sum()} out of {hcopcube.shape[0]}")

    for preg in pregs:
        msk = preg.to_mask()
        slcs_big, slcs_small = msk.get_overlap_slices(mask.shape[1:])

        bmask[insidespecslice, slcs_big[0], slcs_big[1]] |= msk.data.astype('bool')[slcs_small]
        #print(f"total included: {bmask.sum()}")

mhcopcube = hcopcube.with_mask(bmask)
nhcopcube = hcopcube.with_mask(~bmask).with_mask(bmask.any(axis=0))[:,:,250:]
mx_missed = nhcopcube.max(axis=0)
mx_missed.write('hcop_missed_peak_intensity.fits', overwrite=True)

vmax_missed = nhcopcube.with_spectral_unit(u.km/u.s).argmax_world(axis=0)
vmax_missed.write('hcop_missed_peak_velocity.fits', overwrite=True)


import pylab as pl
pl.ion()

fig = pl.figure(1, figsize=(20,20))
fig.clf()
ax1 = pl.subplot(2, 1, 1)

im1 = ax1.imshow(mx_missed.value)
cb1 = pl.colorbar(mappable=im1)
cb1.set_label("K")
ax1.set_title("Peak missed intensity")

ax2 = pl.subplot(2, 1, 2)
im2 = ax2.imshow(vmax_missed.value)
cb2 = pl.colorbar(mappable=im2)
cb2.set_label("km/s")
ax2.set_title("Peak missed velocity")

pl.savefig("peak_missed_hcop.png")
