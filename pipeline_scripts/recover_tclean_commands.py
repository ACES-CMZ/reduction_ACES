"""
Two environmental variables must be set:

    ACES_ROOTDIR: the reduction_ACES directory
    WEBLOG_DIR: the path containing the pipeline* weblogs
"""
import os, sys, glob, json

if os.getenv('ACES_ROOTDIR') is None:
    try:
        os.environ['ACES_ROOTDIR'] = os.path.split(metadata_tools.__file__)[0]
        rootdir = os.environ['ACES_ROOTDIR']
    except ImportError:
        raise ValueError("metadata_tools not found on path; make sure to "
                         "specify ACES_ROOTDIR environment variable "
                         "or your PYTHONPATH variable to include the directory"
                         " containing the ACES code.")
else:
    rootdir = os.environ['ACES_ROOTDIR']
    sys.path.append(rootdir)
    sys.path.append(f'{rootdir}/retrieval_scripts')


from parse_weblog import (get_human_readable_name, get_uid_and_name)


projcode = '2021.1.00172.L'

if 'weblog_dir' not in locals():
    weblog_dir = os.getenv('WEBLOG_DIR')
    if weblog_dir is None:
        raise ValueError("Set an environmental variable 'WEBLOG_DIR' or a local variable `weblog_dir` to specify where to extract the weblogs.")

os.chdir(weblog_dir)


all_cubepars = {}

for pipeline_base in glob.glob(f'{weblog_dir}/pipeline*'):
    sbname, max_baseline = get_human_readable_name(pipeline_base)
    namedict = get_uid_and_name(f"{pipeline_base}/html/t1-1.html")
    mous = namedict['Observing Unit Set Status']
    sbname = namedict['Scheduling Block Name']
    pipeline_run = glob.glob(f"{pipeline_base}/html/stage*")
    print(f"{sbname} = {mous} from {pipeline_base}")
    if 'TP' in sbname:
        print("Skipping TP")
        continue

    cubepars = []
    contpars = []

    cuberun = False
    aggregatecontrun = False
    for stage in pipeline_run:
        logfile = os.path.join(stage, 'casapy.log')
        if os.path.exists(logfile):
            with open(logfile, 'r') as fh:
                for line in fh.readlines():
                    if "hif_makeimlist(specmode='cube')" in line:
                        cuberun=True
                        # skip to next: that will be cubes
                        break
                    if "hif_makeimlist(specmode='cont')" in line:
                        aggrecatecontrun=True
                        # skip to next: that will be the aggregate continuum
                        # (individual continuum is 'mfs')
                        break
                    if 'hif_findcont' in line:
                        # findcont runs tclean but we don't want findcont params
                        break
                    if 'tclean( vis' in line:
                        tcleanpars = line.split("tclean")[-1]
                        tcleanpars = eval(f'dict{tcleanpars}')
                        if tcleanpars['specmode'] == 'cube' and cuberun:
                            cubepars.append(tcleanpars)
                        if tcleanpars['specmode'] == 'mfs':
                            contpars.append(tcleanpars)
            if cubepars and cuberun:
                break
            if contpars and aggregatecontrun:
                break

    if not cubepars or not contpars:
        raise ValueError("No parameters found")

    if not (any('iter1' in pars['imagename'] for pars in cubepars)):
        raise ValueError("TEST")

    # only keep the 'iter1' examples *if* they exist
    cubepars = [pars for pars in cubepars
            if 'iter1' in pars['imagename']]
    contpars = [pars for pars in contpars
            if 'iter1' in pars['imagename']]
    
    # clean out paths from vis
    for pars in cubepars:
        pars["vis"] = [os.path.basename(x) for x in pars["vis"]]
    for pars in contpars:
        pars["vis"] = [os.path.basename(x) for x in pars["vis"]]

    all_cubepars[sbname] = {
            'tclean_cube_pars': cubepars,
            'tclean_cont_pars': contpars,
            'mous': mous,
            'pipeline_run': pipeline_base,
            }

with open(f'{rootdir}/pipeline_scripts/default_tclean_commands.json', 'w') as tcfh:
    json.dump(all_cubepars, tcfh, indent=2)
