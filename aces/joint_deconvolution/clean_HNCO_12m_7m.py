import os
import re
import sys
import glob
import json
import pandas as pd

rootdir = os.getenv('ACES_ROOTDIR')
datadir = os.getenv('ACES_DATA')
workdir = os.getenv('ACES_WORKDIR')
os.chdir(workdir)

with open(f"{rootdir}/aces/pipeline_scripts/default_tclean_commands.json", "r") as fh:
    tclean_commands = json.load(fh)

with open(f"{workdir}/aces_SB_names.csv", "r") as fh:
    sb_names = pd.read_csv(fh)

# Loop over all regions for which we have 12m & 7m data
for i in range(len(sb_names)):
    if not os.path.exists(sb_names['Region'][i]):
        os.mkdir(sb_names['Region'][i])

# Grab all relevant MS files
    visfiles_12m = glob.glob(datadir + '/member.' + sb_names['twelve_m_ID'][i] + '/calibrated/working/*.ms')
    spws_12m     = "25, 27, 29, 31, 33, 35"   # noqa: E221

    visfiles_7m  = glob.glob(datadir + '/member.' + sb_names['seven_m_ID'][i] + '/calibrated/working/*.ms')  # noqa: E221
    spws_7m      = "16, 18, 20, 22, 24, 26"  # noqa: E221

# Split out science target and science SPWs per MS
    if len(visfiles_12m) > 0 and len(visfiles_7m) > 0:
        for visfile in visfiles_12m:
            if not os.path.exists(workdir + '/' + sb_names['Region'][i] + '/' + re.sub(r'^.*?uid___A002', 'uid___A002', visfile).strip('.ms') + '_split.ms'):
                split(vis        = visfile,  # noqa: E221, E251, F821
                      spw        = spws_12m,  # noqa: E221, E251
                      field      = 'Sgr_A_star',  # noqa: E221, E251
                      outputvis  = workdir + '/' + sb_names['Region'][i] + '/' + re.sub(r'^.*?uid___A002', 'uid___A002', visfile).strip('.ms') + '_split.ms',  # noqa: E221, E251
                      datacolumn = 'corrected')  # noqa: E251

        for visfile in visfiles_7m:
            if not os.path.exists(workdir + '/' + sb_names['Region'][i] + '/' + re.sub(r'^.*?uid___A002', 'uid___A002', visfile).strip('.ms') + '_split.ms'):
                split(vis        = visfile,  # noqa: E221, E251, F821
                      spw        = spws_7m,  # noqa: E221, E251
                      field      = 'Sgr_A_star',  # noqa: E221, E251
                      outputvis  = workdir + '/' + sb_names['Region'][i] + '/' + re.sub(r'^.*?uid___A002', 'uid___A002', visfile).strip('.ms') + '_split.ms',  # noqa: E221, E251
                      datacolumn = 'corrected')  # noqa: E251

# Concatenate MSes into a single MS containing only 12 SPWs (6x12m, 6x7m)
        if not os.path.exists(workdir + '/' + sb_names['Region'][i] + '/' + sb_names['Region'][i] + '_12m_7m.concat.ms'):
            concat(vis       = glob.glob(workdir + '/' + sb_names['Region'][i] + '/' + '*split.ms'),  # noqa: E221, E251, F821
                   concatvis = workdir + '/' + sb_names['Region'][i] + '/' + sb_names['Region'][i] + '_12m_7m.concat.ms',  # noqa: E251
                   freqtol   = '10MHz')  # noqa: E221, E251

# Grab 12m cleaning parameters for HNCO SPW > override specific params > clean!
    if (os.path.exists(workdir + '/' + sb_names['Region'][i] + '/' + sb_names['Region'][i] + '_12m_7m.concat.ms') and
            not os.path.exists(workdir + '/' + sb_names['Region'][i] + '/' + sb_names['Region'][i] + '_12m_7m.HNCO.clean.image')):
        hnco_pars                = tclean_commands[sb_names['Reg_original_12m'][i]]['tclean_cube_pars']['spw31']  # noqa: E221
        hnco_pars['vis']         = workdir + '/' + sb_names['Region'][i] + '/' + sb_names['Region'][i] + '_12m_7m.concat.ms'  # noqa: E221
        hnco_pars['spw']         = '3, 9'  # noqa: E221
        hnco_pars['imagename']   = workdir + '/' + sb_names['Region'][i] + '/' + sb_names['Region'][i] + '_12m_7m.HNCO.clean'  # noqa: E221
        hnco_pars['calcpsf']     = True  # noqa: E221
        hnco_pars['calcres']     = True  # noqa: E221
        hnco_pars['cyclefactor'] = 1.5
        hnco_pars['antenna']     = ''  # noqa: E221
        hnco_pars['scan']        = ''  # noqa: E221
        if hnco_pars['nchan'] < 1000:
            hnco_pars['nchan']   = hnco_pars['nchan'] * 2  # noqa: E221
            hnco_pars['width']   = str(float(re.sub('MHz', '', hnco_pars['width'])) / 2.0) + 'MHz'  # noqa: E221

        tclean(**hnco_pars)  # noqa:  F821
