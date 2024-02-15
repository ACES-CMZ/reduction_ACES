import glob
import os
from spectral_cube import SpectralCube
from spectral_cube.lower_dimensional_structures import OneDSpectrum
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy import constants
import pylab as pl
from astropy import stats
from scipy import ndimage

from casatools import ms as mstool

from astropy.io import fits
from astropy.wcs import WCS

from aces.pipeline_scripts.merge_tclean_commands import get_commands
from aces.analysis.parse_contdotdat import parse_contdotdat, cont_channel_selection_to_contdotdat
from aces import conf

import numpy as np
basepath = '/orange/adamginsburg/ACES/'

def make_plot(sbname):


    sbtb = Table.read(f'{basepath}/reduction_ACES/aces/data/tables/aces_SB_uids.csv')
    sbtb.add_index('12m MOUS ID')

    try:
        num = sbtb.loc[sbname].index + 1
    except KeyError:
        print(f"Skipped {sbname} because it's not in the table (maybe it was a reobservation)")
        return


    flagmapfile = fits.open(f'{basepath}/mosaics/continuum/12m_continuum_reimaged_field_number_map.fits')
    flagmap = flagmapfile[0].data[::10, ::10]
    target_wcs = WCS(flagmapfile[0].header)[::10, ::10]

    fig = pl.figure(figsize=(10, 7))
    ax5 = fig.add_subplot(2, 1, 2, projection=target_wcs)
    ax5.contour(flagmap, colors=['k']*int(flagmap.max()), levels=np.arange(flagmap.max()) + 0.5, )

    ax5.contourf(flagmap==num, levels=[0.5, 1.5])

    contdotdatfn = f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibration/cont.dat' 
    if os.path.exists(contdotdatfn):
        contdat = parse_contdotdat(contdotdatfn)
    else:
        print(f"{sbname} has no cont.dat")
        contdat = None

    for ii, spw in enumerate((25, 27, 33, 35)):
        cubefns = (glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibrated/working/*spw{spw}*.statcont.contsub.fits')
                   + glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibrated/working/*sci{spw}*.statcont.contsub.fits')
        )
        if len(cubefns) > 1:
            # filter out s12's
            cubefns = [x for x in cubefns if 's38' in x]
        assert len(cubefns) == 1
        cubefn = cubefns[0]
        cube = SpectralCube.read(cubefn)

        basedir = os.path.dirname(cubefn)
        basename = os.path.splitext(os.path.basename(cubefn))[0]
        specdir = os.path.join(basedir, 'spectra')

        max_fn = f'{specdir}/{basename}.maxspec.fits'
        mean_fn = f'{specdir}/{basename}.meanspec.fits'

        if not os.path.exists(max_fn) or not os.path.exists(mean_fn):
            print(f"Skipping SPW {spw}.  Perhaps some of the max/mean spectra have been run for this field, but not for all windows?")
            continue

        max_spec = fits.getdata(max_fn)
        mean_spec = fits.getdata(mean_fn)

        cont_arr_indiv = np.zeros(cube.shape[0])
        minfreq = cube.spectral_axis.min()
        maxfreq = cube.spectral_axis.max()
        frqarr = u.Quantity(np.linspace(minfreq, maxfreq, cube.shape[0]), u.GHz)
        
        if contdat is not None:
            for frqline in contdat.split(";"):
                fsplit = frqline.split("~")
                f2 = u.Quantity(fsplit[1])
                f1 = u.Quantity(float(fsplit[0]), f2.unit)

                if f1 == f2:
                    continue
                if f1 > f2:
                    f1,f2 = f2,f1

                sel = (frqarr > f1) & (frqarr < f2)
                cont_arr_indiv[sel] = 1

        ax1 = fig.add_subplot(2, 4, ii+1)
        ax1.plot(cube.spectral_axis, max_spec, color='k')
        if contdat is not None:
            max_spec_masked = max_spec.copy()
            max_spec_masked[~cont_arr_indiv.astype('bool')] = np.nan
            ax1.plot(cube.spectral_axis, max_spec_masked, color='r')
        ax1.set_title(str(spw))

        new_contsel = id_continuum(max_spec, threshold)
        max_spec_masked2 = max_spec.copy()
        max_spec_masked2[~new_contsel.astype('bool')] = np.nan
        ax1.plot(cube.spectral_axis, max_spec_masked2, color='lime', alpha=0.75, linewidth=0.25, linestyle=':')

        # ax1.plot(cube.spectral_axis, mean_spec, color='k')
        # mean_spec_masked = mean_spec.copy()
        # mean_spec_masked[~cont_arr_indiv.astype('bool')] = np.nan
        # ax1.plot(cube.spectral_axis, mean_spec_masked, color='r')

    pl.tight_layout()
    pl.savefig(f"{specdir}/{basename}_diagnostic_spectra.png")


def id_continuum(spectrum, threshold=2.5):
    new_contsel = ndimage.binary_dilation(
        ndimage.binary_erosion(
            spectrum < np.nanmedian(spectrum) + 2.5 * stats.mad_std(spectrum),
            iterations=2),
        iterations=1)
    return new_contsel


def assemble_new_contsels():
    cmds = get_commands()

    ms = mstool()

    spwsel = {}

    datadir = f'{conf.basepath}/data/'
    projcode = os.getenv('PROJCODE') or '2021.1.00172.L'
    sous = os.getenv('SOUS') or 'A001_X1590_X30a8'
    gous = os.getenv('GOUS') or 'A001_X1590_X30a9'

    for sbname, allpars in cmds.items():
        if 'TM1' not in sbname:
            # handle 7m separately
            continue
        mous_ = allpars['mous']
        mous = mous_[6:].replace("/", "_")
        assert len(mous) in (14, 15, 16)
        workingpath = f'{datadir}/{projcode}/science_goal.uid___{sous}/group.uid___{gous}/member.uid___{mous}/calibrated/working'
        specdir = f'{workingpath}/spectra/'

        spwsel[sbname] = {'tclean_cont_pars':
                          {'aggregate': {'spw': []},
                           'aggregate_high': {'spw': []}
                            }}

        for msname in allpars['tclean_cont_pars']['aggregate']['vis']:
            vis = f'{workingpath}/{msname}'
            if not os.path.exists(vis):
                vis = vis.replace("targets", "target")
                if not os.path.exists(vis):
                    raise FileNotFoundError(f"{msname} does not exist in {workingpath}")
            ms.open(vis)
            selstrs = []
            selstrs_high = []
            for spw in (25, 27, 33, 35):

                cubefns = (glob.glob(f'{workingpath}/*spw{spw}*.statcont.contsub.fits')
                        + glob.glob(f'{workingpath}/*sci{spw}*.statcont.contsub.fits')
                )
                if len(cubefns) > 1:
                    # filter out s12's
                    cubefns = [x for x in cubefns if 's38' in x]
                assert len(cubefns) == 1
                cubefn = cubefns[0]
                cube = SpectralCube.read(cubefn)

                basename = os.path.splitext(os.path.basename(cubefn))[0]
                max_fn = f'{specdir}/{basename}.maxspec.fits'

                max_spec = OneDSpectrum.from_hdu(fits.open(max_fn))
                contsel_bool = id_continuum(max_spec.value)
                segments, nseg = ndimage.label(contsel_bool)
                frqmins = u.Quantity(ndimage.labeled_comprehension(max_spec.spectral_axis, segments, np.arange(1, nseg+1), np.min, float, np.nan), max_spec.spectral_axis.unit)
                frqmaxs = u.Quantity(ndimage.labeled_comprehension(max_spec.spectral_axis, segments, np.arange(1, nseg+1), np.max, float, np.nan), max_spec.spectral_axis.unit)

                # calculate LSR offset from the data
                native_frqs = ms.cvelfreqs(spw)
                lsr_frqs = ms.cvelfreqs(spw, outframe='LSRK')
                v_offset = np.mean((lsr_frqs - native_frqs) / lsr_frqs) * constants.c

                frqmins_native = (frqmins * (1 - v_offset / constants.c)).to(u.GHz)
                frqmaxs_native = (frqmaxs * (1 - v_offset / constants.c)).to(u.GHz)
                frqpairs = [f'{mn}~{mx}GHz' for mn,mx in zip(frqmins_native.value, frqmaxs_native.value)]

                selstr_ = f'{spw}:' + ";".join(frqpairs)
                selstrs.append(selstr_)
                if spw in (33, 35):
                    selstrs_high.append(selstr_)
            selstr = ",".join(selstrs)
            selstr_high = ",".join(selstrs_high)
            spwsel[sbname]['tclean_cont_pars']['aggregate']['spw'].append(selstr)
            spwsel[sbname]['tclean_cont_pars']['aggregate_high']['spw'].append(selstr_high)

    with open(f"{pipedir}/spw_selections.json", "w") as fh:
        json.dump(fh, spwsel, indent=2)