import numpy as np
from astropy import units as u
from spectral_cube import SpectralCube
from astropy.table import Table
from astropy.io import fits
import regions
from spectral_cube.spectral_cube import _regionlist_to_single_region

# load the HCO+ (and maybe HNCO) cubes
hcopcube = SpectralCube.read('/orange/adamginsburg/cmz/mopra/CMZ_3mm_HCO+.fits')
hncocube = SpectralCube.read('/orange/adamginsburg/cmz/mopra/CMZ_3mm_HNCO.fits')

tbl = Table.read('../SB_naming.tsv', format='ascii.csv', delimiter='\t')

vwidth = 102.585 * 2048 *u.m/u.s

for cube, name in ((hcopcube, 'hcop'), (hncocube, 'hnco')):

    mask = cube.mask.include()
    #print(mask.sum())

    # start with a blank mask and include only observed regions
    bmask = np.zeros_like(mask)


    # create a list of composite regions for labeling purposes
    composites = []

    # loop over the SB_naming table
    for row in tbl:
        indx = row['Proposal ID'][3:]
        # just for debug purposes / progress tracking, print out some stats
        print(f"Row {indx}: mask sum = {bmask.sum()}, which is {bmask.sum() / mask.sum() * 100:0.2f} percent of original")

        # load up regions
        regs = regions.Regions.read(f'../regions/final_cmz{indx}.reg')
        pregs = [reg.to_pixel(cube.wcs.celestial) for reg in regs]

        composite = _regionlist_to_single_region(pregs)
        composite.meta['label'] = indx
        composites.append(composite)

        voff = row['velocity offset']*u.km/u.s
        vmin, vmax = voff-vwidth/2, voff+vwidth/2

        insidespecslice = (cube.spectral_axis > vmin) & (cube.spectral_axis < vmax)
        outsidespecslice = (cube.spectral_axis < vmin) & (cube.spectral_axis > vmax)
        #print(f"Specslice includes {insidespecslice.sum()} and excludes {outsidespecslice.sum()} out of {cube.shape[0]}")

        # loop over all of the individual circle regions in the region file
        for preg in pregs:
            msk = preg.to_mask()
            slcs_big, slcs_small = msk.get_overlap_slices(mask.shape[1:])

            # logic: we're adding by "or" any region that is observed
            bmask[insidespecslice, slcs_big[0], slcs_big[1]] |= msk.data.astype('bool')[slcs_small]
            #print(f"total included: {bmask.sum()}")

    mcube = cube.with_mask(bmask)
    ncube = cube.with_mask(~bmask).with_mask(bmask.any(axis=0))[:,:,250:]
    mx_missed = ncube.max(axis=0)
    mx_missed.write(f'{name}_missed_peak_intensity.fits', overwrite=True)

    vmax_missed = ncube.with_spectral_unit(u.km/u.s).argmax_world(axis=0)
    vmax_missed.write(f'{name}_missed_peak_velocity.fits', overwrite=True)


    import pylab as pl
    pl.ion()

    fig = pl.figure(1, figsize=(14,7))
    fig.clf()
    ax1 = pl.subplot(2, 1, 1, projection=mx_missed.wcs)

    im1 = ax1.imshow(mx_missed.value)
    cb1 = pl.colorbar(mappable=im1)
    cb1.set_label("K")
    ax1.set_title("Peak missed intensity")

    ax2 = pl.subplot(2, 1, 2, projection=vmax_missed.wcs)
    im2 = ax2.imshow(vmax_missed.value)
    cb2 = pl.colorbar(mappable=im2)
    cb2.set_label("km/s")
    ax2.set_title("Peak missed velocity")

    pl.savefig(f"peak_missed_{name}.png", bbox_inches='tight')

    flagmap = np.zeros(cube.shape[1:], dtype='int')

    for comp in composites:
        cmsk = comp.to_mask()

        slcs_big, slcs_small = cmsk.get_overlap_slices(cube.shape[1:])
        flagmap[slcs_big] += (cmsk.data[slcs_small] * int(comp.meta['label'])) * (flagmap[slcs_big] == 0)
        #_, glat, glon = cube.world[0,slcs_big[0],slcs_big[1]]
        #ax2.contour(glon, glat, cmsk.data, levels=[0.5], colors=['r'],
        #           transform=ax2.get_transform('galactic'))

    ax1.contour(flagmap[:,250:], levels=np.arange(flagmap.max()), colors=['w','r']*50,
                linewidths=[0.2]*50)
    ax2.contour(flagmap[:,250:], levels=np.arange(flagmap.max()), colors=['w','r']*50,
                linewidths=[0.2]*50)

    pl.savefig(f"peak_missed_{name}_regcontours.png", bbox_inches='tight')


    fig2 = pl.figure(2)
    pl.imshow(flagmap[:,250:])

    flagmapfits = fits.PrimaryHDU(data=flagmap[:,250:],
                                  header=cube.wcs.celestial[:,250:].to_header())
    flagmapfits.writeto('region_number_map.fits', overwrite=True)



    voff = 0
    vmin, vmax = voff-vwidth/2, voff+vwidth/2

    outsidespecslice = (cube.spectral_axis < vmin) | (cube.spectral_axis > vmax)

    origslab = cube.with_mask(bmask.any(axis=0)).with_mask(outsidespecslice[:,None,None])[:,:,250:]
    mxmissed_noshift = origslab.max(axis=0)
    vmaxmissed_noshift = origslab.with_spectral_unit(u.km/u.s).argmax_world(axis=0)

    fig = pl.figure(3, figsize=(14,7))
    fig.clf()
    ax1 = pl.subplot(2, 1, 1, projection=mxmissed_noshift.wcs)

    im1 = ax1.imshow(mxmissed_noshift.value)
    cb1 = pl.colorbar(mappable=im1)
    cb1.set_label("K")
    ax1.set_title("Peak missed intensity")

    ax2 = pl.subplot(2, 1, 2, projection=vmaxmissed_noshift.wcs)
    im2 = ax2.imshow(vmaxmissed_noshift.value)
    cb2 = pl.colorbar(mappable=im2)
    cb2.set_label("km/s")
    ax2.set_title("Peak missed velocity")

    pl.savefig(f"peak_missed_{name}_noshift.png", bbox_inches='tight')

    ax1.contour(flagmap[:,250:], levels=np.arange(flagmap.max()), colors=['w','r']*50,
                linewidths=[0.2]*50)
    ax2.contour(flagmap[:,250:], levels=np.arange(flagmap.max()), colors=['w','r']*50,
                linewidths=[0.2]*50)

    pl.savefig(f"peak_missed_{name}_noshift_regcontours.png", bbox_inches='tight')
