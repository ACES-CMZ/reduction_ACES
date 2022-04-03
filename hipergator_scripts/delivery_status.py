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

def get_mous_to_sb_mapping(project_code):
    from astroquery.alma import Alma

    tbl = Alma.query(payload={'project_code': project_code}, cache=False,
                     public=False)['member_ous_uid','schedblock_name', 'qa2_passed']
    mapping = {row['member_ous_uid']: row['schedblock_name'] for row in tbl if row['qa2_passed'] == 'T'}
    return mapping

mousmap = get_mous_to_sb_mapping('2021.1.00172.L')
mousmap_ = {key.replace("/","_").replace(":","_"):val for key,val in mousmap.items()}

cwd = os.getcwd()
basepath = '/orange/adamginsburg/ACES/rawdata/'
os.chdir(basepath)

datatable = {}

spwlist = {'12M': {33, 35, 25, 27, 29, 31},
           'TM1': {33, 35, 25, 27, 29, 31},
           '7M': {16, 18, 20, 22, 24, 26}}

def wildexists(x):
    return len(glob.glob(x)) > 0

basepath = dataroot = '/orange/adamginsburg/ACES/data/2021.1.00172.L'
#workpath = '/blue/adamginsburg/adamginsburg/ACES/workdir/'

for fullpath in glob.glob(f"{basepath}/sci*/group*/member*/"):
    mous = os.path.basename(fullpath.strip('/')).split(".")[-1]
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

            for spw in spwlist[config]:
                # /orange/adamginsburg/ACES/rawdata/2021.1.00172.L/
                # science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_Xae/calibrated/working/uid___A001_X15a0_Xae.s9_0.Sgr_A_star_sci.spw26.mfs.I.iter1.image
                bn = f'{mous}.s*_0.Sgr_A_star_sci.spw{spw}.{clean}.I'
                workingpath = f'{fullpath}/calibrated/working/'

                imageglob = f'{workingpath}/{bn}.image'
                pbcorglob = f'{workingpath}/{bn}.image.pbcor'
                psfglob = f'{workingpath}/{bn}.psf'

                exists = (wildexists(pbcorglob)
                          or ("WIPim" if wildexists(imageglob)
                              else "WIPpsf" if wildexists(psfglob)
                              else False))


                if field not in datatable:
                    datatable[field] = {}
                if spw not in datatable[field]:
                    datatable[field][spw] = {}
                if clean not in datatable[field][spw]:
                    datatable[field][spw][clean] = {}


                datatable[field][spw][clean] = {
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


os.chdir(cwd)
