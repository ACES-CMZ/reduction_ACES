import numpy as np
import itertools
import warnings
import glob
import os
from astropy.io import fits
from astropy import visualization
from astropy.table import Table, Column
from spectral_cube import SpectralCube
from astropy.stats import mad_std
from astropy import log
import pylab as pl
import subprocess
from aces.retrieval_scripts.mous_map import get_mous_to_sb_mapping
from aces import conf

def get_mousmap_(**kwargs):
    mousmap = get_mous_to_sb_mapping('2021.1.00172.L', **kwargs)
    mousmap_ = {key.replace("/","_").replace(":","_"):val for key,val in mousmap.items()}
    return mousmap_


cwd = os.getcwd()
basepath = conf.basepath
os.chdir(basepath)

datatable = {}

spwlist = {'12M': {25, 27, 29, 31, 33, 35, },
           'TM1': {25, 27, 29, 31, 33, 35, },
           '7M': {16, 18, 20, 22, 24, 26},
           'TP':  {16, 18, 20, 22, 24, 26},
           }

def wildexists(x):
    return len(glob.glob(x)) > 0

def main():
    datapath = dataroot = f'{basepath}/data/2021.1.00172.L'
    #workpath = '/blue/adamginsburg/adamginsburg/ACES/workdir/'

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
            raise ValueError(f"sbname={sbname} is not being handled correctly by delivery_status.py")

        if ' ' in config:
            # handle this case: 'Sgr_A_st_aj_03_7M Sgr_A_st_aj_03_7M_original'
            config = config.split(" ")[0]
        rerun = 'original' in sbname

        if 'updated' in sbname:
            field = field+"_updated"
        if 'original' in sbname:
            field = field+"_original"

        for clean in ('mfs', 'cont', 'cube',):
            for suffix in (".image",  ):#".contsub.image"):#, ".contsub.JvM.image.fits", ".JvM.image.fits"):

                #globblob = f'{fullpath}/calibrated/working/*{clean}*.iter1{suffix}'
                #fn = glob.glob(f'{dataroot}/{globblob}')

                for spwn in sorted(spwlist[config] | ({'aggregate'} if clean!='cube' else set()), key=lambda x: str(x)):
                    # /orange/adamginsburg/ACES/rawdata/2021.1.00172.L/
                    # science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_Xae/calibrated/working/uid___A001_X15a0_Xae.s9_0.Sgr_A_star_sci.spw26.mfs.I.iter1.image
                    spw = spwn if isinstance(spwn, str) else f'spw{spwn}'

                    # we use 'spwkey' because the glob search string is different from the adopted name
                    # (we search for spw_*.cont, but the name is 'aggregate)
                    spwkey = spw

                    # aggregate continuum is named with the full list of spws
                    if spw == 'aggregate':
                        spw = 'spw*'
                        clean_ = 'cont' # _not_ mfs
                    else:
                        clean_ = clean

                    bn = f'{mous}.s*_0.Sgr_A_star_sci.{spw}.{clean_}.I'
                    workingpath = f'{fullpath}/calibrated/working/'


                    tts = '.tt0' if spwkey == 'aggregate' else ''

                    imageglob = f'{workingpath}/{bn}.iter1.image{tts}'
                    pbcorglob = f'{workingpath}/{bn}.iter1.image{tts}.pbcor'
                    psfglob = f'{workingpath}/{bn}.iter1.psf{tts}'

                    exists = (wildexists(pbcorglob)
                              or ("WIPim" if wildexists(imageglob)
                                  else "WIPpsf" if wildexists(psfglob)
                                  else False))


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
    for sbname,xx in datatable.items():
        cols = {}
        tables[sbname] = cols
        for config, yy in xx.items():
            rows = []
            for spw, zz in yy.items():
                for imtype, ww in zz.items():
                    rows = ww
                    #cols[f'{config}.{spw}'].append(','.join([argtype for argtype, status in ww.items() if status is True]))
                cols[f'{config}.{spw}'] = rows
        rows = [[cn,] + list(v.values()) for cn,v in cols.items()]
        tables[sbname] = Table(rows=rows, names=['config.spw'] + list(ww.keys()))

    with open(f'{basepath}/reduction_ACES/aces/data/tables/imaging_completeness_grid.txt', 'w') as fh:
        for key,tb in tables.items():
            fh.write(f'SB {key}:\n')
            fh.write("\n".join(tb.pformat()))
            fh.write("\n\n")


    os.chdir(cwd)

    # Hack for debugging
    globals().update(locals())
    return locals()
