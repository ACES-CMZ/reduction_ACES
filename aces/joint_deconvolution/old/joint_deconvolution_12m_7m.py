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

with open(f"{rootdir}/aces/joint_deconvolution/aces_SB_names.csv", "r") as fh:
    sb_names = pd.read_csv(fh)

with open(f"{rootdir}/aces/pipeline_scripts/default_tclean_commands.json", "r") as fh:
    tclean_commands = json.load(fh)

hnco_pars = tclean_commands[region + '_TM1']['tclean_cube_pars']['spw31']

if not os.path.isdir(pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.image') and not os.path.isdir(pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_TP.SPW3'):
    tclean(vis=pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3_concat.ms',
           field='Sgr_A_star',
           datacolumn='corrected',
           imagename=region + '_12m_7m_SPW3.clean',
           specmode='cube',
           gridder='mosaic',
           restoringbeam='common',
           calcpsf=True,
           calcres=True,
           restoration=True,
           pbcor=True,
           parallel=True,
           interactive=0,
           cyclefactor=1.5,
           mosweight=hnco_pars['mosweight'],
           outframe=hnco_pars['outframe'],
           intent=hnco_pars['intent'],
           imsize=hnco_pars['imsize'],
           cell=hnco_pars['cell'],
           phasecenter=hnco_pars['phasecenter'],
           perchanweightdensity=hnco_pars['perchanweightdensity'],
           width=hnco_pars['width'],
           start=hnco_pars['start'],
           nchan=hnco_pars['nchan'],
           deconvolver=hnco_pars['deconvolver'],
           weighting=hnco_pars['weighting'],
           robust=hnco_pars['robust'],
           niter=hnco_pars['niter'],
           threshold=hnco_pars['threshold'],
           usemask=hnco_pars['usemask'],
           sidelobethreshold=hnco_pars['sidelobethreshold'],
           noisethreshold=hnco_pars['noisethreshold'],
           lownoisethreshold=hnco_pars['lownoisethreshold'],
           negativethreshold=hnco_pars['negativethreshold'],
           minbeamfrac=hnco_pars['minbeamfrac'],
           growiterations=hnco_pars['growiterations'],
           dogrowprune=hnco_pars['dogrowprune'])
