import glob
import os
from astropy import log
from aces.retrieval_scripts.mous_map import get_mous_to_sb_mapping
from aces import conf
import time


def get_mousmap_(**kwargs):
    mousmap = get_mous_to_sb_mapping('2021.1.00172.L', **kwargs)
    mousmap_ = {key.replace("/", "_").replace(":", "_"): val
                for key, val in mousmap.items()}
    return mousmap_


def wildexists(x):
    return len(glob.glob(x)) > 0


def main():
    t0 = time.time()
    cwd = os.getcwd()
    basepath = conf.basepath
    os.chdir(basepath)

    datatable = {}

    spwlist = {'12M': [25, 27, 29, 31, 33, 35],
               'TM1': [25, 27, 29, 31, 33, 35],
               '7M': [16, 18, 20, 22, 24, 26],
               'TP': [16, 18, 20, 22, 24, 26],
               }

    datapath = dataroot = f'{basepath}/data/2021.1.00172.L'
    # workpath = '/blue/adamginsburg/adamginsburg/ACES/workdir/'

    looplist = glob.glob(f"{datapath}/sci*/group*/member*/")
    assert len(looplist) > 0

    looplist = sorted(looplist, key=lambda x: os.path.basename(x))

    mousmap_ = get_mousmap_()

    assert len(mousmap_) > 0

    for fullpath in looplist:
        log.info(f'Working on {fullpath}')
        mous = os.path.basename(fullpath.strip('/')).split(".")[-1]
        try:
            sbname = mousmap_[mous]
        except KeyError:
            mousmap_ = get_mousmap_(refresh=True)
            sbname = mousmap_[mous]
        field = sbname.split("_")[3]

        # we have a few cases:
        # 'Sgr_A_st_j_03_TM1_updated'
        # 'Sgr_A_st_b_updated_03_7M'
        if sbname.split("_")[4] == 'updated' or sbname.split("_")[5] == 'updated':
            config = sbname.split("_")[6]
        else:
            config = sbname.split("_")[5]

        # sanity check
        if config in ('updated', 'original', '03'):
            raise ValueError(
                f"sbname={sbname} is not being handled correctly by delivery_status.py")

        if ' ' in config:
            # handle this case: 'Sgr_A_st_aj_03_7M Sgr_A_st_aj_03_7M_original'
            config = config.split(" ")[0]
        rerun = 'original' in sbname

        if 'updated' in sbname:
            field = field + "_updated"
        if 'original' in sbname:
            field = field + "_original"

        # 'cont' mode isn't really used, even though it's more proper
        # 'mfs' is remapped to 'cont' in the job runner and in the renaming
        # if-statement below
        for clean in ('mfs', 'cube',):  # 'cont',
            for suffix in (".image",):
                # ".contsub.image"):#, ".contsub.JvM.image.fits", ".JvM.image.fits"):

                # globblob = f'{fullpath}/calibrated/working/*{clean}*.iter1{suffix}'
                # fn = glob.glob(f'{dataroot}/{globblob}')

                if clean == 'cube':
                    spwlist_ = spwlist[config]
                else:
                    spwlist_ = ['aggregate', 'aggregate_high', 'aggregate_low']

                for spwn in sorted(spwlist_, key=lambda x: str(x)):
                    # /orange/adamginsburg/ACES/rawdata/2021.1.00172.L/
                    # science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_Xae/calibrated/working/uid___A001_X15a0_Xae.s9_0.Sgr_A_star_sci.spw26.mfs.I.iter1.image
                    spw = spwn if isinstance(spwn, str) else f'spw{spwn}'

                    # we use 'spwkey' because the glob search string is different from the adopted name
                    # (we search for spw_*.cont, but the name is 'aggregate)
                    spwkey = spw

                    # the name in the files is cont, not mfs, for aggregate
                    clean_ = 'cont'  # _not_ mfs
                    # aggregate continuum is named with the full list of spws
                    if spw == 'aggregate':
                        spw = "spw" + "_".join(str(x) for x in spwlist[config])
                    elif spw == 'aggregate_high':
                        spw = "spw" + "_".join(str(x) for x in spwlist[config][-2:])
                    elif spw == 'aggregate_low':
                        spw = "spw" + "_".join(str(x) for x in spwlist[config][:2])
                    else:
                        clean_ = clean

                    for stepglob in ("s*_0.", ""):
                        for scidot in ('sci.', 'sci'):
                            for spwtxt in ('spw', ''):
                                bn = f'{mous}.{stepglob}Sgr_A_star_{scidot}{spw}.{clean_}.I'.replace('spw', spwtxt)
                                workingpath = f'{fullpath}/calibrated/working/'

                                tts = '.tt0' if 'aggregate' in spwkey else ''

                                for iter_or_manual in ('iter1', 'manual', 'iter1.reclean', 'manual.reclean', 'iter2.reclean'):
                                    imageglob = f'{workingpath}/{bn}.{iter_or_manual}.image{tts}'
                                    pbcorglob = f'{workingpath}/{bn}.{iter_or_manual}.image{tts}.pbcor'
                                    psfglob = f'{workingpath}/{bn}.{iter_or_manual}.psf{tts}'

                                    exists = (wildexists(pbcorglob) or
                                              ("WIPim" if wildexists(imageglob)
                                              else "WIPpsf" if wildexists(psfglob)
                                              else False))
                                    if exists:
                                        break
                                if exists:
                                    break
                            if exists:
                                break
                        if exists:
                            break
                    if 'X184' in mous:
                        print(mous, bn, exists)

                    if mous not in datatable:
                        datatable[mous] = {}
                    if config not in datatable[mous]:
                        datatable[mous][config] = {}
                    if spwkey not in datatable[mous][config]:
                        datatable[mous][config][spwkey] = {}
                    if clean not in datatable[mous][config][spwkey]:
                        datatable[mous][config][spwkey][clean] = {}

                    datatable[mous][config][spwkey][clean] = {
                        "field": field,
                        "image": wildexists(imageglob),
                        "pbcor": wildexists(pbcorglob),
                        "psf": wildexists(psfglob),
                        "WIP": ("Done" if wildexists(pbcorglob)
                                else "WIPim" if wildexists(imageglob)
                                else "WIPpsf" if wildexists(psfglob)
                                else False)
                    }

    import json
    with open(f'{basepath}/reduction_ACES/aces/data/tables/imaging_completeness_grid.json', 'w') as fh:
        json.dump(datatable, fh)

    # make a table
    from astropy.table import Table
    tables = {}
    for sbname, xx in datatable.items():
        cols = {}
        tables[sbname] = cols
        for config, yy in xx.items():
            rows = []
            for spw, zz in yy.items():
                for imtype, ww in zz.items():
                    rows = ww
                    # cols[f'{config}.{spw}'].append(','.join([argtype for argtype, status in ww.items() if status is True]))
                cols[f'{config}.{spw}'] = rows
        rows = [[cn, ] + list(v.values()) for cn, v in cols.items()]
        tables[sbname] = Table(
            rows=rows, names=['config.spw'] + list(ww.keys()))

    with open(f'{basepath}/reduction_ACES/aces/data/tables/imaging_completeness_grid.txt', 'w') as fh:
        for key, tb in tables.items():
            fh.write(f'SB {key}:\n')
            fh.write("\n".join(tb.pformat()))
            fh.write("\n\n")

    os.chdir(cwd)

    t1 = time.time()
    print(f"delivery_status took {t1 - t0} seconds = {(t1 - t0) / 60} minutes = {(t1 - t0) / 3600} hours")

    # Hack for debugging
    globals().update(locals())
    return locals()
