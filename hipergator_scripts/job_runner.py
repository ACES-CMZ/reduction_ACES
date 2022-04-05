import numpy as np
import os
import json
import glob
import shutil
import copy
from mous_map import get_mous_to_sb_mapping

basepath = "/orange/adamginsburg/ACES/rawdata"
grouppath = f"{basepath}/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9"

mouses = [os.path.basename(x)
        for x in
        glob.glob(f'{grouppath}/member.uid___A001_X15a0_X*')]

mousmap = get_mous_to_sb_mapping('2021.1.00172.L')
mousmap_ = {key.replace("/","_").replace(":","_"):val for key,val in mousmap.items()}

parameters = {'member.uid___A001_X15a0_Xea': 
        {'mem': 128, 'ntasks': 32, 'mpi': True,  } }
newpars = parameters.copy()

default_parameters = {f'{os.path.basename(mous.strip("/"))}':
                      {'mem': 128, 'ntasks': 32, 'mpi': True, }
    for mous in mouses}

for key in newpars:
    default_parameters[key] = newpars[key]

parameters = default_parameters

if __name__ == "__main__":
    import subprocess
    import datetime
    import os
    import json
    from astropy.io import ascii
    import sys

    verbose = '--verbose' in sys.argv

    with open('/orange/adamginsburg/web/secure/ACES/tables/imaging_completeness_grid.json', 'r') as fh:
        imaging_status = json.load(fh)

    sacct = subprocess.check_output(['sacct',
                                   '--format=JobID,JobName%45,Account%15,QOS%17,State,Priority%8,ReqMem%8,CPUTime%15,Elapsed%15,Timelimit%15,NodeList%20'])
    tbl = ascii.read(sacct.decode().split("\n"))

    scriptpath = '/orange/adamginsburg/ACES/reduction_ACES/hipergator_scripts/'

    qos = os.getenv('QOS')
    if qos:
        account = os.environ['ACCOUNT'] = 'adamginsburg' if 'adamginsburg' in qos else 'astronomy-dept'
        if 'astronomy-dept' not in qos and 'adamginsburg' not in qos:
            raise ValueError(f"Invalid QOS {qos}")
    else:
        account = 'astronomy-dept'
        qos = 'astronomy-dept-b'
    logpath = os.environ['LOGPATH']='/blue/adamginsburg/adamginsburg/ACES/logs/'


    for mous,spwpars in parameters.items():

        sbname = mousmap_[mous.split('.')[1]]
        field = sbname.split("_")[3]
        config = sbname.split("_")[5]

        do_contsub = bool(spwpars.get('do_contsub'))
        contsub_suffix = '.contsub' if do_contsub else ''

        for config in imaging_status[field]:
            for spw in imaging_status[field][config]:
                for imtype in imaging_status[field][config][spw]:
                    imstatus = imaging_status[field][config][spw][imtype]
                    if imstatus['image'] and imstatus['pbcor']:
                        # Done
                        continue
                    elif imstatus['WIP']:
                        if '--redo-wip' in sys.argv:
                            print(f"field {field} {spw} {imtype} is in progress: imstatus={imstatus['WIP']}; trying anyway (if it is not in the 'PENDING' or 'RUNNING' queue)")
                        else:
                            continue
                    else:
                        if verbose:
                            print(f"field {field} {spw} {imtype} has not begun: imstatus={imstatus}")

                    if verbose:
                        print(f"spw={spw}, spwpars={spwpars}")


                    workdir = '/blue/adamginsburg/adamginsburg/ACES/workdir'
                    jobname = f"{field}_{config}_{spw}_{imtype}"

                    match = tbl['JobName'] == jobname
                    if any(match):
                        states = np.unique(tbl['State'][match])
                        if 'RUNNING' in states:
                            jobid = tbl['JobID'][match & (tbl['State'] == 'RUNNING')]
                            continue
                            print(f"Skipped job {jobname} because it's RUNNING as {set(jobid)}")
                        elif 'PENDING' in states:
                            jobid = tbl['JobID'][match & (tbl['State'] == 'PENDING')]
                            print(f"Skipped job {jobname} because it's PENDING as {set(jobid)}")
                            continue
                        elif 'COMPLETED' in states:
                            jobid = tbl['JobID'][match & (tbl['State'] == 'COMPLETED')]
                            if '--redo-completed' in sys.argv:
                                print(f"Redoing job {jobname} even though it's COMPLETED as {set(jobid)} (if it is not pending)")
                            else:
                                print(f"Skipped job {jobname} because it's COMPLETED as {set(jobid)}")
                                continue
                        elif 'FAILED' in states:
                            jobid = tbl['JobID'][match & (tbl['State'] == 'FAILED')]
                            if '--redo-failed' in sys.argv:
                                print(f"Redoing job {jobname} even though it's FAILED as {set(jobid)}")
                            else:
                                print(f"Skipped job {jobname} because it's FAILED as {set(jobid)}")
                                continue
                        elif 'TIMEOUT' in states:
                            jobid = tbl['JobID'][match & (tbl['State'] == 'TIMEOUT')]
                            print(f"Restarting job {jobname} because it TIMED OUT as {set(jobid)}")


                    # handle specific parameters
                    mem = int(spwpars["mem"])
                    os.environ['MEM'] = mem = f'{mem}gb'
                    ntasks = spwpars["ntasks"]
                    os.environ['NTASKS'] = str(ntasks)
                    os.environ['DO_CONTSUB'] = str(do_contsub)
                    os.environ['SLURM_NTASKS'] = str(ntasks)
                    os.environ['WORK_DIRECTORY'] = workdir
                    os.environ['FIELD_ID'] = field

                    #/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/
                    #group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X130/calibrated/working/
                    #tclean_cube_pars_Sgr_A_st_y_03_TM1_spw27.py

                    calwork = f'{grouppath}/{mous}/calibrated/working'
                    cleantype = {'cube': 'cube', 'mfs': 'cont'}[imtype]
                    try:
                        spwname = f'spw{int(spw)}'
                    except Exception:
                        spwname = spw
                    scriptname = f'{calwork}/tclean_{cleantype}_pars_Sgr_A_st_{field}_03_{config}_{spwname}.py'
                    if (imtype == 'mfs' and spw == 'aggregate') or imtype == 'cube':
                        assert os.path.exists(scriptname)
                    else:
                        # skip MFS individual spws
                        continue
                    os.environ['SCRIPTNAME'] = scriptname

                    basename = f'{field}_{spw}_{imtype}{contsub_suffix}'
                    # basename = "{0}_{1}_spw{2}_{3}".format(field, band, spw, arrayname)

                    # it is safe to remove things beyond here because at this point we're committed
                    # to re-running
                    if '--dry-run' not in sys.argv:
                        if '--remove-failed' in sys.argv:
                            #print(f"Removing files matching '{workdir}/{basename}.*'")
                            failed_files = glob.glob(f'{workdir}/{basename}.*')
                            if any('.image' in x for x in failed_files):
                                print(f"Found a .image in the failed file list: {failed_files}.  Continuing.")
                                #raise ValueError(f"Found a .image in the failed file list: {failed_files}")
                            else:
                                for ff in failed_files:
                                    print(f"Removing {ff}")
                                    shutil.rmtree(ff)
                        
                        tempdir_name = f'{field}_{spw}_{imtype}_{contsub_suffix}'
                        print(f"Removing files matching '{workdir}/{tempdir_name}/IMAGING_WEIGHT.*'")
                        old_tempfiles = (glob.glob(f'{workdir}/{tempdir_name}/IMAGING_WEIGHT*') +
                                         glob.glob(f'{workdir}/{tempdir_name}/TempLattice*'))
                        for tfn in old_tempfiles:
                            print(f"Removing {tfn}")
                            shutil.rmtree(tfn)

                    if spwpars['mpi']:
                        mpisuffix = '_mpi'
                        cpus_per_task = 1
                        os.environ['SLURM_TASKS_PER_NODE'] = str(ntasks)
                    else:
                        assert ntasks == 1
                        mpisuffix = ''
                        cpus_per_task = ntasks

                    os.environ['CPUS_PER_TASK'] = str(cpus_per_task)

                    runcmd = f'{scriptpath}/run_simple_script.sh'

                    now = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
                    os.environ['LOGFILENAME'] = f"{logpath}/casa_log_line_{jobname}_{now}.log"


                    cmd = f'sbatch --ntasks={ntasks} --cpus-per-task={cpus_per_task} --mem={mem} --output={jobname}_%j.log --job-name={jobname} --account={account} --qos={qos} --export=ALL  {runcmd}'

                    if '--dry-run' in sys.argv:
                        if verbose:
                            print(cmd)
                        print()
                        #print(subprocess.check_output('env').decode())
                        #raise
                    else:
                        sbatch = subprocess.check_output(cmd.split())

                        print(f"Started sbatch job with jobid={sbatch.decode()} and parameters {spwpars} and script {scriptname}")
