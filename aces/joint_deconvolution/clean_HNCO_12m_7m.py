# flake8: noqa
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

with open(f"{rootdir}/aces/joint_deconvolution/aces_SB_names.csv", "r") as fh:
    sb_names = pd.read_csv(fh)

for i in range(len(sb_names)):
    if not os.path.exists(sb_names['Region name'][i]):
        os.mkdir(sb_names['Region name'][i])

    visfiles_12m = glob.glob(datadir+'/member.uid___'+sb_names['twelve_m'][i]+'/calibrated/working/*.ms')
    spws_12m     = "25, 27, 29, 31, 33, 35"

    visfiles_7m  = glob.glob(datadir+'/member.uid___'+sb_names['seven_m'][i]+'/calibrated/working/*.ms')
    spws_7m      = "16, 18, 20, 22, 24, 26"

    if len(visfiles_12m) > 0 and len(visfiles_7m) > 0:
        for visfile in visfiles_12m:
            if not os.path.exists(workdir+'/'+sb_names['Region name'][i]+'/'+re.sub(r'^.*?uid___A002', 'uid___A002', visfile).strip('.ms')+'_split.ms'):
                split(vis        = visfile,
                      spw        = spws_12m,
                      field      = 'Sgr_A_star',
                      outputvis  = workdir+'/'+sb_names['Region name'][i]+'/'+re.sub(r'^.*?uid___A002', 'uid___A002', visfile).strip('.ms')+'_split.ms',
                      datacolumn ='corrected')

        for visfile in visfiles_7m:
            if not os.path.exists(workdir+'/'+sb_names['Region name'][i]+'/'+re.sub(r'^.*?uid___A002', 'uid___A002', visfile).strip('.ms')+'_split.ms'):
                split(vis        = visfile,
                      spw        = spws_7m,
                      field      = 'Sgr_A_star',
                      outputvis  = workdir+'/'+sb_names['Region name'][i]+'/'+re.sub(r'^.*?uid___A002', 'uid___A002', visfile).strip('.ms')+'_split.ms',
                      datacolumn ='corrected')

        if not os.path.exists(workdir+'/'+sb_names['Region name'][i]+'/'+sb_names['Region name'][i]+'_12m_7m.concat.ms'):
            concat(vis       = glob.glob(workdir+'/'+sb_names['Region name'][i]+'/'+'*split.ms'),
                   concatvis = workdir+'/'+sb_names['Region name'][i]+'/'+sb_names['Region name'][i]+'_12m_7m.concat.ms',
                   freqtol   = '10MHz')

    if os.path.exists(workdir+'/'+sb_names['Region name'][i]+'/'+sb_names['Region name'][i]+'_12m_7m.concat.ms') and not os.path.exists(workdir+'/'+sb_names['Region name'][i]+'/'+sb_names['Region name'][i]+'_12m_7m.HNCO.clean'):
        hnco_pars = tclean_commands[sb_names['Region name'][i]+'_03_TM1']['tclean_cube_pars']['spw31']
        tclean(vis                  = workdir+'/'+sb_names['Region name'][i]+'/'+sb_names['Region name'][i]+'_12m_7m.concat.ms',
               field                = 'Sgr_A_star',
               datacolumn           = 'corrected',
               spw                  = '3, 9',
               imagename            = workdir+'/'+sb_names['Region name'][i]+'/'+sb_names['Region name'][i]+'_12m_7m.HNCO.clean',
               specmode             = 'cube',
               gridder              = 'mosaic',
               restoringbeam        = 'common',
               calcpsf              = True,
               calcres              = True,
               restoration          = True,
               pbcor                = True,
               parallel             = True,
               interactive          = 0,
               cyclefactor          = 2.0,
               mosweight            = hnco_pars['mosweight'],
               outframe             = hnco_pars['outframe'],
               intent               = hnco_pars['intent'],
               imsize               = hnco_pars['imsize'],
               cell                 = hnco_pars['cell'],
               phasecenter          = hnco_pars['phasecenter'],
               perchanweightdensity = hnco_pars['perchanweightdensity'],
               width                = hnco_pars['width'],
               start                = hnco_pars['start'],
               nchan                = hnco_pars['nchan'],
               deconvolver          = hnco_pars['deconvolver'],
               weighting            = hnco_pars['weighting'],
               robust               = hnco_pars['robust'],
               niter                = hnco_pars['niter'],
               threshold            = hnco_pars['threshold'],
               usemask              = hnco_pars['usemask'],
               sidelobethreshold    = hnco_pars['sidelobethreshold'],
               noisethreshold       = hnco_pars['noisethreshold'],
               lownoisethreshold    = hnco_pars['lownoisethreshold'],
               negativethreshold    = hnco_pars['negativethreshold'],
               minbeamfrac          = hnco_pars['minbeamfrac'],
               growiterations       = hnco_pars['growiterations'],
               dogrowprune          = hnco_pars['dogrowprune'])
