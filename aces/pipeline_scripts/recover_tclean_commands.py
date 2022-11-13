"""
Two environmental variables must be set:

    ACES_ROOTDIR: the reduction_ACES directory
    WEBLOG_DIR: the path containing the pipeline* weblogs
"""
from aces.retrieval_scripts.parse_weblog import (get_human_readable_name, get_uid_and_name)
import re
import os
import glob
import json
from astropy import log

from aces import conf
rootdir = conf.basepath

# if os.getenv('ACES_ROOTDIR') is None:
#     try:
#         os.environ['ACES_ROOTDIR'] = os.path.split(metadata_tools.__file__)[0]
#         rootdir = os.environ['ACES_ROOTDIR']
#     except ImportError:
#         raise ValueError("metadata_tools not found on path; make sure to "
#                          "specify ACES_ROOTDIR environment variable "
#                          "or your PYTHONPATH variable to include the directory"
#                          " containing the ACES code.")
# else:
#     rootdir = os.environ['ACES_ROOTDIR']
#     sys.path.append(rootdir)
#     sys.path.append(f'{rootdir}/retrieval_scripts')


def main():

    projcode = '2021.1.00172.L'

    source_name = 'Sgr_A_star'

    if 'weblog_dir' not in locals():
        weblog_dir = os.getenv('WEBLOG_DIR')
        if weblog_dir is None:
            raise ValueError("Set an environmental variable 'WEBLOG_DIR' or a local variable `weblog_dir` to specify where to extract the weblogs.")

    os.chdir(weblog_dir)

    cube_re = re.compile(r"hif_makeimlist\(.*specmode='cube'")
    mfs_re = re.compile(r"hif_makeimlist\(.*specmode='mfs'")
    cont_re = re.compile(r"hif_makeimlist\(.*specmode='cont'")

    all_cubepars = {}

    for pipeline_base in glob.glob(f'{weblog_dir}/pipeline*'):
        sbname, max_baseline = get_human_readable_name(pipeline_base, verbose=False)
        # DEBUG tool - uncomment as needed
        # if '_al_' not in sbname:
        #     print(".", end="")
        #     continue
        # print()
        namedict = get_uid_and_name(f"{pipeline_base}/html/t1-1.html")
        mous = namedict['Observing Unit Set Status']
        sbname = namedict['Scheduling Block Name']
        pipeline_run = glob.glob(f"{pipeline_base}/html/stage*")

        if len(pipeline_run) < 26:
            log.info(f"Skipping pipeline {pipeline_base} from mous {mous} sb {sbname}"
                     f"because it only had {len(pipeline_run)} stages and was therefore likely a rerun.")
            continue

        # need them to run in order
        pipeline_run = sorted(pipeline_run, key=lambda x: int(x.split("stage")[-1]))

        if 'TP' in sbname:
            log.info(f"Skipping TP {sbname} = {mous} from {pipeline_base}:  TP Skipped!")
            continue
        else:
            log.info(f"Processing {sbname} = {mous} from {pipeline_base}")

        cubepars = {}
        contpars = {}

        cuberun = False
        aggregatecontrun = False
        for stage in pipeline_run:
            logfile = os.path.join(stage, 'casapy.log')
            if os.path.exists(logfile):
                with open(logfile, 'r') as fh:
                    for line in fh.readlines():
                        if cube_re.search(line):
                            # can be
                            # hif_makeimlist(specmode='cube')
                            # or
                            # hif_makeimlist(specmode='cont', robust=1.0)
                            cuberun = True
                            # skip to next: that will be cubes
                            break
                        if cont_re.search(line):
                            aggregatecontrun = True
                            # skip to next: that will be the aggregate continuum
                            # (individual continuum is 'mfs')
                            break
                        if 'hif_findcont' in line:
                            # findcont runs tclean but we don't want findcont params
                            break
                        if 'tclean( vis' in line:
                            log.debug(f"{logfile} found tclean command")
                            tcleanpars = line.split("tclean")[-1]
                            tcleanpars = eval(f'dict{tcleanpars}')
                            if tcleanpars['field'] == source_name:
                                if tcleanpars['specmode'] == 'cube' and cuberun:
                                    spw = set(tcleanpars['spw'])
                                    if len(spw) != 1:
                                        raise ValueError("Cube found multiple spws")
                                    else:
                                        # get the string value out
                                        spw = next(iter(spw))
                                    cubepars[f'spw{spw}'] = tcleanpars
                                if tcleanpars['specmode'] == 'mfs' and aggregatecontrun:
                                    contpars['aggregate'] = tcleanpars
                log.debug(f"{logfile}: aggregate continuum = {aggregatecontrun}, cube = {cuberun}")
                if cubepars and cuberun:
                    cuberun = False
                if contpars and aggregatecontrun:
                    aggregatecontrun = False

        if not cubepars or not contpars:
            if sbname in ('Sgr_A_st_al_03_TM1',):
                print(f"Skipping {sbname} because it was manually QA'd")
                continue
            raise ValueError("No parameters found")

        if not (any('iter1' in pars['imagename'] for pars in cubepars.values())):
            print("There is no iter1 in the imagename parameters")
            print(f"sbname {sbname} has images {[pars['imagename'] for pars in cubepars.values()]}")
            print("Skipping")
            continue

        # only keep the 'iter1' examples *if* they exist
        cubepars = {key: pars for key, pars in cubepars.items()
                    if 'iter1' in pars['imagename']}
        contpars = {key: pars for key, pars in contpars.items()
                    if 'iter1' in pars['imagename']}

        # clean out paths from vis
        for pars in cubepars.values():
            pars["vis"] = [os.path.basename(x) for x in pars["vis"]]
        for pars in contpars.values():
            pars["vis"] = [os.path.basename(x) for x in pars["vis"]]

        log.debug(f'For sb {sbname} = {mous}, added parameters')
        all_cubepars[sbname] = {
            'tclean_cube_pars': cubepars,
            'tclean_cont_pars': contpars,
            'mous': mous,
            'pipeline_run': pipeline_base,
        }

        # DEBUG tool.  Uncomment as needed
        # if '_al_' in sbname and 'TM1' in sbname:
        #     globals().update(locals())
        #     return all_cubepars

    with open(f'{rootdir}/reduction_ACES/aces/pipeline_scripts/default_tclean_commands.json', 'w') as tcfh:
        json.dump(all_cubepars, tcfh, indent=2)

    globals().update(locals())

if __name__ == "__main__":
    main()
