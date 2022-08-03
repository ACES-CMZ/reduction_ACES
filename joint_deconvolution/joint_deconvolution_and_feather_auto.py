import json, os, sys, glob

region      = sys.argv[1]
pth         = os.path.abspath(os.getcwd())

with open(f"./default_tclean_commands.json", "r") as fh:
    tclean_commands = json.load(fh)

os.chdir(pth+'/'+region)
if not os.path.exists('clean_12m_7m'):
    os.mkdir('clean_12m_7m')
os.chdir(pth+'/'+region+'/clean_12m_7m/')

vis_all     = glob.glob(pth+'/'+region+'/twelve_m/spws/*.spw3') + glob.glob(pth+'/'+region+'/seven_m/spws/*.spw3')
if os.path.isdir(region+'_12m_7m_SPW3_concat.ms')==False:
    concat(vis       = vis_all,
           concatvis = region+'_12m_7m_SPW3_concat.ms')

hnco_pars   = tclean_commands[region+'_TM1']['tclean_cube_pars']['spw31']

if os.path.isdir(pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.image')==False and os.path.isdir(pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_TP.SPW3')==False:
    tclean(vis                     = pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3_concat.ms',
           field                   = 'Sgr_A_star',
           datacolumn              = 'corrected',
           imagename               = region+'_12m_7m_SPW3.clean',
           specmode                = 'cube',
           gridder                 = 'mosaic',
           restoringbeam           = 'common',
           calcpsf                 = True,
           calcres                 = True,
           restoration             = True,
           pbcor                   = True,
           parallel                = True,
           interactive             = 0,
           cyclefactor             = 1.5,
           mosweight               = hnco_pars['mosweight'],
           outframe                = hnco_pars['outframe'],
           intent                  = hnco_pars['intent'],
           imsize                  = hnco_pars['imsize'],
           cell                    = hnco_pars['cell'],
           phasecenter             = hnco_pars['phasecenter'],
           perchanweightdensity    = hnco_pars['perchanweightdensity'],
           width                   = hnco_pars['width'],
           start                   = hnco_pars['start'],
           nchan                   = hnco_pars['nchan'],
           deconvolver             = hnco_pars['deconvolver'],
           weighting               = hnco_pars['weighting'],
           robust                  = hnco_pars['robust'],
           niter                   = hnco_pars['niter'],
           threshold               = hnco_pars['threshold'],
           usemask                 = hnco_pars['usemask'],
           sidelobethreshold       = hnco_pars['sidelobethreshold'],
           noisethreshold          = hnco_pars['noisethreshold'],
           lownoisethreshold       = hnco_pars['lownoisethreshold'],
           negativethreshold       = hnco_pars['negativethreshold'],
           minbeamfrac             = hnco_pars['minbeamfrac'],
           growiterations          = hnco_pars['growiterations'],
           dogrowprune             = hnco_pars['dogrowprune'])

tp_data = glob.glob(pth+'/'+region+'/total_power/'+'*.spw23.cube.I.sd.fits')[0]

if os.path.isdir(pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.image')==True:
    tp_freq       = imhead(tp_data,mode='get',hdkey='restfreq')
    high_res_freq = imhead(pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.image',mode='get',hdkey='restfreq')

    if tp_freq['value'] != high_res_freq['value']:
        imreframe(imagename = tp_data,
                  restfreq  = high_res_freq['value']+' Hz',
                  output    = tp_data+'.reframe')

    if os.path.isdir(tp_data+'.reframe')==True:
        imregrid(imagename  = tp_data+'.reframe',
                 template   = pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.image',
                 output     = tp_data+'.regrid')
    else:
        imregrid(imagename  = tp_data,
                 template   = pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.image',
                 output     = tp_data+'.regrid')

    imsubimage(imagename    = pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.image',
               outfile      = pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.image.dropdeg',
               dropdeg      = True)

    imsubimage(imagename    = pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.pb',
               outfile      = pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.pb.dropdeg',
               dropdeg      = True)

    imsubimage(imagename    = tp_data+'.regrid',
               outfile      = tp_data+'.regrid.dropdeg',
               dropdeg      = True)

    immath(imagename        = [tp_data+'.regrid.dropdeg',
                               pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.pb.dropdeg'],
           expr             = 'IM0*IM1',
           outfile          = tp_data+'.regrid.dropdeg.depb')

    feather(imagename       = pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_TP.SPW3',
            highres         = pth+'/'+region+'/clean_12m_7m/'+region+'_12m_7m_SPW3.clean.image.dropdeg',
            lowres          = tp_data+'.regrid.dropdeg.depb')
