import os
import sys
import glob
import string
import numpy as np
from casatools import quanta
sys.path.append('scripts')
from parse_contdotdat import parse_contdotdat, contchannels_to_linechannels

ms = mstool()
qq = quanta()

region = sys.argv[1]
pth    = os.path.abspath(os.getcwd())

os.chdir(pth+'/'+region+'/twelve_m/')
if not os.path.exists('spws'):
    os.mkdir('spws')

visfiles = glob.glob('*.ms')
contfile = 'cont.dat'
spws     = [25, 27, 29, 31, 33, 35]

freqs = {}
ms.open(visfiles[0])
for spw in spws:
    freqs[spw] = ms.cvelfreqs(spwids=[spw], outframe='LSRK')
ms.close()

cont_channel_selection = parse_contdotdat(contfile)
linechannels, linefracs = contchannels_to_linechannels(cont_channel_selection,
                                                       freqs,
                                                       return_fractions=True)

for visfile in visfiles:
    if os.path.isdir(visfile.strip('.ms')+'_target.ms')==False:
        split(vis        = visfile,
              field      = 'Sgr_A_star',
              intent     = "OBSERVE_TARGET#ON_SOURCE",
              outputvis  = visfile.strip('.ms')+'_target.ms',
              datacolumn = 'corrected')

    if os.path.isdir(visfile.strip('.ms')+'_target.ms.contsub')==False:
        uvcontsub(vis           = visfile.strip('.ms')+'_target.ms',
                  fitspw        = linechannels,
                  field         = 'Sgr_A_star',
                  excludechans  = True,
                  combine       = 'spw',
                  solint        = 'int',
                  fitorder      = 0,
                  want_cont     = False)

    i=0
    for spw in spws:
        if os.path.isdir('./spws/'+visfile.strip('.ms')+'.spw'+str(i))==False:
            split(vis        = visfile.strip('.ms')+'_target.ms.contsub',
                  spw        = spw,
                  field      = 'Sgr_A_star',
                  outputvis  = './spws/'+visfile.strip('.ms')+'.spw'+str(i),
                  datacolumn = 'data')
        i=i+1

os.chdir(pth)

