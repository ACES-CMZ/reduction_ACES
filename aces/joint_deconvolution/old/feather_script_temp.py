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

tp_data = glob.glob(pth + '/' + region + '/total_power/' + '*.spw23.cube.I.sd.fits')[0]

if os.path.isdir(pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.image'):
    tp_freq = imhead(tp_data, mode='get', hdkey='restfreq')
    high_res_freq = imhead(pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.image', mode='get', hdkey='restfreq')

    if tp_freq['value'] != high_res_freq['value']:
        imreframe(imagename=tp_data,
                  restfreq=high_res_freq['value'] + ' Hz',
                  output=tp_data + '.reframe')

    if os.path.isdir(tp_data + '.reframe'):
        imregrid(imagename=tp_data + '.reframe',
                 template=pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.image',
                 output=tp_data + '.regrid')
    else:
        imregrid(imagename=tp_data,
                 template=pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.image',
                 output=tp_data + '.regrid')

    imsubimage(imagename=pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.image',
               outfile=pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.image.dropdeg',
               dropdeg=True)

    imsubimage(imagename=pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.pb',
               outfile=pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.pb.dropdeg',
               dropdeg=True)

    imsubimage(imagename=tp_data + '.regrid',
               outfile=tp_data + '.regrid.dropdeg',
               dropdeg=True)

    immath(imagename=[tp_data + '.regrid.dropdeg',
                      pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.pb.dropdeg'],
           expr='IM0*IM1',
           outfile=tp_data + '.regrid.dropdeg.depb')

    feather(imagename=pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_TP.SPW3',
            highres=pth + '/' + region + '/clean_12m_7m/' + region + '_12m_7m_SPW3.clean.image.dropdeg',
            lowres=tp_data + '.regrid.dropdeg.depb')
