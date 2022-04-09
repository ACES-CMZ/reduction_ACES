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
from mous_map import get_mous_to_sb_mapping

def get_mousmap_(**kwargs):
    mousmap = get_mous_to_sb_mapping('2021.1.00172.L', **kwargs)
    mousmap_ = {key.replace("/","_").replace(":","_"):val for key,val in mousmap.items()}
    return mousmap_

mousmap_ = get_mousmap_()

cwd = os.getcwd()
basepath = '/orange/adamginsburg/ACES/rawdata/'
os.chdir(basepath)

datatable = {}

spwlist = {'12M': {25, 27, 29, 31, 33, 35, },
           'TM1': {25, 27, 29, 31, 33, 35, },
           '7M': {16, 18, 20, 22, 24, 26},
           'TP':  {16, 18, 20, 22, 24, 26},
           }

def wildexists(x):
    return len(glob.glob(x)) > 0

basepath = dataroot = '/orange/adamginsburg/ACES/data/2021.1.00172.L'
#workpath = '/blue/adamginsburg/adamginsburg/ACES/workdir/'

looplist = glob.glob(f"{basepath}/sci*/group*/member*/")
looplist = sorted(looplist, key=lambda x: os.path.basename(x))

for fullpath in looplist:
    mous = os.path.basename(fullpath.strip('/')).split(".")[-1]
    try:
        sbname = mousmap_[mous]
    except KeyError:
        mousmap_ = get_mousmap_(refresh=True)
        sbname = mousmap_[mous]
    field = sbname.split("_")[3]
    config = sbname.split("_")[5]
    if ' ' in config:
        # handle this case: 'Sgr_A_st_aj_03_7M Sgr_A_st_aj_03_7M_original'
        config = config.split(" ")[0]
    rerun = 'original' in sbname

    for clean in ('mfs', 'cube'):
        for suffix in (".image", ):#".contsub.image"):#, ".contsub.JvM.image.fits", ".JvM.image.fits"):
            #globblob = f'{fullpath}/calibrated/working/*{clean}*.iter1{suffix}'
            #fn = glob.glob(f'{dataroot}/{globblob}')

            for spwn in sorted(spwlist[config] | {'aggregate'} if clean=='mfs' else set(), key=lambda x: str(x)):
                # /orange/adamginsburg/ACES/rawdata/2021.1.00172.L/
                # science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_Xae/calibrated/working/uid___A001_X15a0_Xae.s9_0.Sgr_A_star_sci.spw26.mfs.I.iter1.image
                spw = spwn if isinstance(spwn, str) else f'spw{spwn}'
                bn = f'{mous}.s*_0.Sgr_A_star_sci.{spw}.{clean}.I'
                workingpath = f'{fullpath}/calibrated/working/'

                imageglob = f'{workingpath}/{bn}.iter1.image'
                pbcorglob = f'{workingpath}/{bn}.iter1.image.pbcor'
                psfglob = f'{workingpath}/{bn}.iter1.psf'

                exists = (wildexists(pbcorglob)
                          or ("WIPim" if wildexists(imageglob)
                              else "WIPpsf" if wildexists(psfglob)
                              else False))


                if field not in datatable:
                    datatable[field] = {}
                if config not in datatable[field]:
                    datatable[field][config] = {}
                if spw not in datatable[field][config]:
                    datatable[field][config][spw] = {}
                if clean not in datatable[field][config][spw]:
                    datatable[field][config][spw][clean] = {}


                datatable[field][config][spw][clean] = {
                        "image": wildexists(imageglob),
                        "pbcor": wildexists(pbcorglob),
                        "psf": wildexists(psfglob),
                        "WIP": ("Done" if wildexists(pbcorglob)
                                else "WIPim" if wildexists(imageglob)
                                else "WIPpsf" if wildexists(psfglob)
                                else False)
                        }

import json
with open('/orange/adamginsburg/web/secure/ACES/tables/imaging_completeness_grid.json', 'w') as fh:
    json.dump(datatable, fh)

# make a table
from astropy.table import Table
tables = {}
for field,xx in datatable.items():
    cols = {}
    tables[field] = cols
    for config, yy in xx.items():
        rows = []
        for spw, zz in yy.items():
            for imtype, ww in zz.items():
                rows = ww
                #cols[f'{config}.{spw}'].append(','.join([argtype for argtype, status in ww.items() if status is True]))
            cols[f'{config}.{spw}'] = rows
    rows = [[cn,] + list(v.values()) for cn,v in cols.items()]
    tables[field] = Table(rows=rows, names=['config.spw'] + list(ww.keys()))

with open('/orange/adamginsburg/web/secure/ACES/tables/imaging_completeness_grid.txt', 'w') as fh:
    for key,tb in tables.items():
        fh.write(f'Field {key}:\n')
        fh.write("\n".join(tb.pformat()))
        fh.write("\n\n")


os.chdir(cwd)
