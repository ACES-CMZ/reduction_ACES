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

from aces.pipeline_scripts.merge_tclean_commands import get_commands
from aces.analysis.parse_contdotdat import parse_contdotdat, cont_channel_selection_to_contdotdat
from aces import conf

import json

import numpy as np
basepath = '/orange/adamginsburg/ACES/'
pipedir = os.path.realpath(os.path.dirname(__file__) + "/../pipeline_scripts/")


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
    ax5.contour(flagmap, colors=['k'] * int(flagmap.max()), levels=np.arange(flagmap.max()) + 0.5, )

    ax5.contourf(flagmap == num, levels=[0.5, 1.5])

    contdotdatfn = f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibration/cont.dat'
    if os.path.exists(contdotdatfn):
        contdat = parse_contdotdat(contdotdatfn)
    else:
        print(f"{sbname} has no cont.dat")
        contdat = None

    for ii, spw in enumerate((25, 27, 33, 35)):
        cubefns = (glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibrated/working/*spw{spw}*.cube.I.iter1.image')
                   + glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibrated/working/*sci{spw}*.cube.I.iter1.image')
                   + glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibrated/working/*sci{spw}*.cube.I.manual.image')
                   + glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibrated/working/*spw{spw}*.cube.I.manual.image')
                   + glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibrated/working/*sci{spw}*.cube.I.iter1.reclean.image')
                   + glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibrated/working/*spw{spw}*.cube.I.manual.reclean.image')
                   + glob.glob(f'{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_{sbname}/calibrated/working/*spw{spw}*.cube.I.iter1.reclean.image')
        )
        if len(cubefns) > 1:
            # filter out s12's
            cubefns = [x for x in cubefns if 's38' in x]
        if len(cubefns) == 0:
            print(f"NO CUBE FOUND FOR {sbname} {spw}")
            print(f"NO CUBE FOUND FOR {sbname} {spw}")
            print(f"NO CUBE FOUND FOR {sbname} {spw}")
            print(f"NO CUBE FOUND FOR {sbname} {spw}")
            print(f"NO CUBE FOUND FOR {sbname} {spw}")
            print(f"NO CUBE FOUND FOR {sbname} {spw}")
            continue
        assert len(cubefns) == 1
        cubefn = cubefns[0]
        cube = SpectralCube.read(cubefn)

        basedir = os.path.dirname(cubefn)
        # splitext is needed if we use .fits files
        # basename = os.path.splitext(os.path.basename(cubefn))[0]
        basename = os.path.basename(cubefn)
        if cubefn.endswith('.image'):
            assert basename.endswith('.image')
        specdir = os.path.join(basedir, 'spectra')

        max_fn = f'{specdir}/{basename}.maxspec.fits'
        mean_fn = f'{specdir}/{basename}.meanspec.fits'

        if not os.path.exists(max_fn) or not os.path.exists(mean_fn):
            print(max_fn, os.path.exists(max_fn))
            print(mean_fn, os.path.exists(mean_fn))
            print(f"Skipping SPW {spw}.  Perhaps some of the max/mean spectra have been run for this field, but not for all windows?")
            #raise
            continue

        max_spec = fits.getdata(max_fn)
        # mean_spec = fits.getdata(mean_fn)

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
                    f1, f2 = f2, f1

                sel = (frqarr > f1) & (frqarr < f2)
                cont_arr_indiv[sel] = 1

        ax1 = fig.add_subplot(2, 4, ii + 1)
        ax1.plot(cube.spectral_axis, max_spec, color='k')
        if contdat is not None:
            max_spec_masked = max_spec.copy()
            max_spec_masked[~cont_arr_indiv.astype('bool')] = np.nan
            ax1.plot(cube.spectral_axis, max_spec_masked, color='r')

        new_contsel = id_continuum(max_spec, threshold=2.5)
        max_spec_masked2 = max_spec.copy()
        max_spec_masked2[~new_contsel.astype('bool')] = np.nan
        ax1.plot(cube.spectral_axis, max_spec_masked2, color='lime',
                 alpha=0.75, linewidth=0.75, linestyle=':')

        # include the number of channels and total bandwidth in the plot titles
        tot_bw = np.mean(np.diff(frqarr)) * new_contsel.sum()
        ax1.set_title(f'{spw}: {new_contsel.sum()}ch {tot_bw:0.1f}')
        # zoom the plot in to show the selection
        mn = np.nanmin(max_spec_masked2)
        std = np.nanstd(max_spec_masked2)
        ax1.set_ylim(mn - 5 * std,
                     np.nanmax(max_spec_masked2) + 5 * std)
        if new_contsel.sum() == 0:
            raise ValueError("Found no continuum")

        # ax1.plot(cube.spectral_axis, mean_spec, color='k')
        # mean_spec_masked = mean_spec.copy()
        # mean_spec_masked[~cont_arr_indiv.astype('bool')] = np.nan
        # ax1.plot(cube.spectral_axis, mean_spec_masked, color='r')

    pl.tight_layout()
    pl.savefig(f"{specdir}/{basename}_diagnostic_spectra.png")
    pl.close(fig.number)


def id_continuum(spectrum, threshold=2.5):
    med = np.nanmedian(spectrum)
    mad = stats.mad_std(spectrum, ignore_nan=True)

    new_contsel = ndimage.binary_dilation(
        ndimage.binary_erosion(
            ((spectrum < med + threshold * mad) &
             (spectrum > med - threshold * mad)),
            iterations=2),
        iterations=1)
    return new_contsel


def assemble_new_contsels(convert_to_native=False, allow_missing_maxspec=False):
    cmds = get_commands()

    ms = mstool()

    spwsel = {}

    datadir = f'{conf.basepath}/data/'
    projcode = os.getenv('PROJCODE') or '2021.1.00172.L'
    sous = os.getenv('SOUS') or 'A001_X1590_X30a8'
    gous = os.getenv('GOUS') or 'A001_X1590_X30a9'

    for sbname, allpars in cmds.items():
        mous_ = allpars['mous']
        mous = mous_[6:].replace("/", "_")
        assert len(mous) in (14, 15, 16)
        workingpath = f'{datadir}/{projcode}/science_goal.uid___{sous}/group.uid___{gous}/member.uid___{mous}/calibrated/working'
        specdir = f'{workingpath}/spectra/'

        spwsel[sbname] = {'tclean_cont_pars':
                          {'aggregate': {'spw': []},
                           'aggregate_high': {'spw': []},
                           'aggregate_low': {'spw': []}
                            }}

        for msname in allpars['tclean_cont_pars']['aggregate']['vis']:
            vis = f'{workingpath}/{msname}'
            if not os.path.exists(vis):
                vis = vis.replace("targets", "target")
                if not os.path.exists(vis):
                    vis = vis.replace("_targets", "")
                    vis = vis.replace("_target", "")
                    if not os.path.exists(vis):
                        raise FileNotFoundError(f"{msname} does not exist in {workingpath}")
            ms.open(vis)
            selstrs = []
            selstrs_high = []
            selstrs_low = []

            all_cubes = glob.glob(f'{workingpath}/*cube*image')
            first_spw = int(all_cubes[0].split("spw")[-1].split("sci")[-1][:2])
            if first_spw in (16, 18, 20, 22, 24, 26):
                array = '7m'
                spwset = (16, 18, 24, 26)
            elif first_spw in (25, 27, 29, 31, 33, 35):
                array = '12m'
                spwset = (25, 27, 33, 35)
            else:
                print(f"first_spw = {first_spw}, which is not a known SPW")
                continue

            for spw in spwset:

                cubefns = [x if f'sci{spw}' in x or f'spw{spw}' in x
                           for x in all_cubes]

                if len(cubefns) > 1:
                    # filter out s12's
                    cubefns = [x for x in cubefns if 's38' in x]
                if len(cubefns) == 0:
                    print(f"NO CUBE FOUND FOR {sbname} {spw}")
                    print(f"NO CUBE FOUND FOR {sbname} {spw}")
                    print(f"NO CUBE FOUND FOR {sbname} {spw}")
                    print(f"NO CUBE FOUND FOR {sbname} {spw}")
                    print(f"NO CUBE FOUND FOR {sbname} {spw}")
                    print(f"NO CUBE FOUND FOR {sbname} {spw}")
                    continue
                assert len(cubefns) == 1
                cubefn = cubefns[0]
                # cube = SpectralCube.read(cubefn)

                if cubefn.endswith('.image'):
                    basename = os.path.basename(cubefn)
                elif cubefn.endswith('.fits'):
                    basename = os.path.splitext(os.path.basename(cubefn))[0]
                else:
                    raise ValueError("Unrecognized file type")
                max_fn = f'{specdir}/{basename}.maxspec.fits'
                if not os.path.exists(max_fn):
                    print(f"NO MAX SPECTRUM FOUND FOR {sbname} {spw}")
                    print(f"NO MAX SPECTRUM FOUND FOR {sbname} {spw}")
                    print(f"NO MAX SPECTRUM FOUND FOR {sbname} {spw}")
                    print(f"NO MAX SPECTRUM FOUND FOR {sbname} {spw}")
                    print(f"NO MAX SPECTRUM FOUND FOR {sbname} {spw}")
                    print(f"NO MAX SPECTRUM FOUND FOR {sbname} {spw}")
                    if allow_missing_maxspec:
                        continue
                    else:
                        raise ValueError("No max spectrum found.")

                max_spec = OneDSpectrum.from_hdu(fits.open(max_fn))
                contsel_bool = id_continuum(max_spec.value)
                segments, nseg = ndimage.label(contsel_bool)
                frqmins = u.Quantity(ndimage.labeled_comprehension(max_spec.spectral_axis, segments, np.arange(1, nseg + 1), np.min, float, np.nan), max_spec.spectral_axis.unit)
                frqmaxs = u.Quantity(ndimage.labeled_comprehension(max_spec.spectral_axis, segments, np.arange(1, nseg + 1), np.max, float, np.nan), max_spec.spectral_axis.unit)

                if convert_to_native:
                    # calculate LSR offset from the data
                    native_frqs = ms.cvelfreqs(spw)
                    lsr_frqs = ms.cvelfreqs(spw, outframe='LSRK')
                    v_offset = np.mean((lsr_frqs - native_frqs) / lsr_frqs) * constants.c

                    frqmins_native = (frqmins * (1 - v_offset / constants.c)).to(u.GHz)
                    frqmaxs_native = (frqmaxs * (1 - v_offset / constants.c)).to(u.GHz)
                    frqpairs = [f'{mn}~{mx}GHz' for mn, mx in zip(frqmins_native.value, frqmaxs_native.value)]
                else:
                    frqpairs = [f'{mn}~{mx}GHz' for mn, mx in zip(frqmins.value, frqmaxs.value)]

                selstr_ = f'{spw}:' + ";".join(frqpairs)
                selstrs.append(selstr_)
                if spw in (33, 35):
                    selstrs_high.append(selstr_)
                elif spw in (25, 27):
                    selstrs_low.append(selstr_)
            selstr = ",".join(selstrs)
            selstr_high = ",".join(selstrs_high)
            selstr_low = ",".join(selstrs_low)
            spwsel[sbname]['tclean_cont_pars']['aggregate']['spw'].append(selstr)
            spwsel[sbname]['tclean_cont_pars']['aggregate_high']['spw'].append(selstr_high)
            spwsel[sbname]['tclean_cont_pars']['aggregate_low']['spw'].append(selstr_low)

    with open(f"{pipedir}/spw_selections.json", "w") as fh:
        json.dump(spwsel, fh, indent=2)
