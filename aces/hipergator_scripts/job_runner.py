import numpy as np
import os
import json
import glob
import shutil
import subprocess
import datetime
from astropy.io import ascii
import sys
from astropy import log
from aces.retrieval_scripts.mous_map import get_mous_to_sb_mapping
from aces import conf

basepath = conf.basepath
datapath = f"{basepath}/data"
grouppath = f"{datapath}/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9"

mouses = [os.path.basename(x)
          for x in
          glob.glob(f'{grouppath}/member.uid___A001_X15*_X*')]

parameters = {'member.uid___A001_X15a0_Xea': {'mem': 128, 'ntasks': 32, 'mpi': True, },  # the usual MPI crash error is occurring
              'member.uid___A001_X15a0_X142': {'mem': 128, 'ntasks': 1, 'mpi': False, },  # ditto
              'member.uid___A001_X15a0_Xca': {'mem': 128, 'ntasks': 1, 'mpi': False, },  # field ag: MPI crash
              'member.uid___A001_X15a0_X160': {'mem': 128, 'ntasks': 1, 'mpi': False, },
              'member.uid___A001_X15a0_X1a2': {'mem': 256, 'ntasks': 1, 'mpi': False, },  # field ar: timeout
              'member.uid___A001_X15a0_X1a8': {'mem': 256, 'ntasks': 1, 'mpi': False, },  # write lock frozen spw33 OOMs; try MPI?.  MPI write-locks everything.
              'member.uid___A001_X15a0_Xa6': {'mem': 256, 'ntasks': 64, 'mpi': True, },  # spw33 is taking for-ev-er; try MPI?  created backup first: backup_20221108_beforempi/
              'member.uid___A001_X15a0_X190': {'mem': 256, 'ntasks': 1, 'mpi': False,
                                               'jobtime': '200:00:00', 'burst': False},  # ao: same as above, too long.  But, MPI fails with writelock. NON-MPI also fails!?
              'member.uid___A001_X15a0_X14e': {'mem': 256, 'ntasks': 64, 'mpi': True, },  # ad: same as above, too long
              'member.uid___A001_X15a0_Xd0': {'mem': 256, 'ntasks': 1, 'mpi': False, },  # field i spw35: timeout
              }
newpars = parameters.copy()

# June 1, 2022: try using fewer tasks to see if it reduces likelihood of race condition
# Idea based on CASA log: try using nprocs/2 + 1 MPI services
# July 14, 2022: the failure rate is ~1, so let's just say 'f it'
default_parameters = {f'{os.path.basename(mous.strip("/"))}':
                      {'mem': 128, 'ntasks': 1, 'mpi': False, }
                      for mous in mouses}

for key in newpars:
    default_parameters[key] = newpars[key]

parameters = default_parameters


def main():

    verbose = '--verbose' in sys.argv
    debug = '--debug' in sys.argv
    check_syntax = '--check-syntax' in sys.argv

    if debug:
        log.setLevel('DEBUG')

    with open(f'{basepath}/reduction_ACES/aces/data/tables/imaging_completeness_grid.json', 'r') as fh:
        imaging_status = json.load(fh)

    mousmap = get_mous_to_sb_mapping('2021.1.00172.L')
    mousmap_ = {key.replace("/", "_").replace(":", "_"): val for key, val in mousmap.items()}

    sacct = subprocess.check_output(['/opt/slurm/bin/sacct',
                                     '--format=JobID,JobName%45,Account%15,QOS%17,State,Priority%8,ReqMem%8,CPUTime%15,Elapsed%15,Timelimit%15,NodeList%20'])
    tbl = ascii.read(sacct.decode().split("\n"))

    scriptpath = f'{basepath}/reduction_ACES/aces/hipergator_scripts/'

    qos = os.getenv('QOS')
    if qos:
        account = os.environ['ACCOUNT'] = 'adamginsburg' if 'adamginsburg' in qos else 'astronomy-dept'
        if 'astronomy-dept' not in qos and 'adamginsburg' not in qos:
            raise ValueError(f"Invalid QOS {qos}")
    else:
        account = 'astronomy-dept'
        qos = 'astronomy-dept-b'
    logpath = os.environ['LOGPATH'] = conf.logpath

    for mous, spwpars in parameters.items():
        mousname = mous.split('.')[1]

        jobtime = '96:00:00'
        if 'burst' in spwpars:
            if 'jobtime' in spwpars:
                jobtime = spwpars.pop('jobtime')
            if not spwpars.pop('burst'):
                qos_ = qos.strip('-b')
                print(f"Using non-burst QOS {qos_} with job time {jobtime}")
        else:
            qos_ = qos

        if mousname not in mousmap_:
            log.error(f"mousname {mousname}  (from {mous}) is not in the mousmap.")
            continue

        sbname = mousmap_[mousname]
        field = sbname.split("_")[3]
        config = sbname.replace("_updated", "").split("_")[5]

        assert config in ('7M', 'TM1', 'TP')

        do_contsub = bool(spwpars.get('do_contsub'))
        contsub_suffix = '.contsub' if do_contsub else ''

        log.debug(f"mous={mous} field={field} sbname={sbname} config={config}")

        if mousname not in imaging_status:
            log.warn(f"WARNING: MOUS {mousname} aka {mous} (field={field}, sbname={sbname}, config={config}) has no delivery status; it may be incomplete")
            continue

        for config_ in imaging_status[mousname]:
            log.debug(f"mous={mous} field={field} sbname={sbname} config={config} config_={config_}")

            if config_ == 'TP':
                # we don't do TP
                log.debug(f"TOTAL POWER, SKIPPING: field={field}, config={config}, sbname={sbname}, mous={mousname}")
                continue
            if config_ != config:
                # imaging_status doesn't know which config is being asked for
                # skip if the config is not the right one for the mous
                log.debug(f"Loop iteration config {config_} is not the requested config {config}, SKIPPING")
                continue
            if not os.path.exists(f'{grouppath}/{mous}'):
                print(f"MOUS {mousname} is not downloaded/extracted (path={grouppath}/{mous}).")
                continue
            for spw in imaging_status[mousname][config]:

                for imtype in imaging_status[mousname][config][spw]:
                    log.debug(f"spw={spw} imtype={imtype}{'**************AGGREGATE**********' if 'aggregate' in imtype else ''}")
                    imstatus = imaging_status[mousname][config][spw][imtype]

                    calwork = f'{grouppath}/{mous}/calibrated/working'
                    cleantype = {'cube': 'cube', 'mfs': 'cont'}[imtype]
                    try:
                        spwname = f'spw{int(spw)}'
                    except Exception:
                        spwname = spw
                    scriptname = f'{calwork}/tclean_{cleantype}_pars_Sgr_A_st_{field}_03_{config}_{spwname}.py'
                    scriptname_glob = f'{calwork}/tclean_{cleantype}_pars_Sgr_A_st_{field}*_03_{config}*_{spwname}.py'
                    if (imtype == 'mfs' and 'aggregate' in spw) or imtype == 'cube':
                        if not os.path.exists(scriptname):
                            if any(glob.glob(scriptname_glob)):
                                scriptname = glob.glob(scriptname_glob)[0]
                            else:
                                if '7M' in scriptname and imtype != 'cube':
                                    # April 22, 2023: I don't think I ever made aggregate high scripts
                                    # for the 7m data?  So the error message is erroneous / not useful
                                    log.debug(f"Skipped {scriptname} b/c it's 7M aggregate")
                                    continue
                                else:
                                    print(f"ERROR: script {scriptname} does not exist!  This may indicate that `write_tclean_scripts` has not been run or the pipeline hasn't finished.")
                                if mousname in 'member.uid___A001_X15b4_X3d':
                                    raise ValueError("The script definitely _does_ exist.")
                                continue
                        else:
                            # all is good
                            pass
                    else:
                        # skip MFS individual spws
                        log.debug(f"imtype is {imtype} and spw is {spw}.  SKIPPING because it's an MFS single-window")
                        continue

                    if check_syntax:
                        # check syntax
                        flake = subprocess.run(['/orange/adamginsburg/miniconda3/envs/python39/bin/flake8',
                                                '--select=E999',
                                                scriptname],
                                               check=True)
                        if flake.returncode == 0:
                            log.debug(flake.stdout)
                        else:
                            raise SyntaxError(flake.stdout)

                    os.environ['SCRIPTNAME'] = scriptname

                    if imstatus['image'] and imstatus['pbcor']:
                        # Done
                        log.debug(f"imtype={imtype} spw={spw}: Image status is done! (image & pbcor).  SKIPPING")
                        continue
                    elif imstatus['WIP']:
                        if '--redo-wip' in sys.argv:
                            print(f"field {field} {spw} {imtype} is in progress: imstatus={imstatus['WIP']}; trying anyway (if it is not in the 'PENDING' or 'RUNNING' queue)")
                        else:
                            log.debug(f"Image status is WIP:{imstatus['WIP']}.  --redo-wip was not specified.   SKIPPING")
                            continue
                    else:
                        if verbose:
                            print(f"field {field} {spw} {imtype} has not begun: imstatus={imstatus}")

                    if verbose:
                        print(f"spw={spw}, imtype={imtype}, spwpars={spwpars}, imstatus={imstatus}")

                    workdir = conf.workpath
                    # TODO: configure job_runner & write_tclean_scripts to work with SLURM_TMPDIR as the working directory?
                    # workdir = "SLURM_TMPDIR"
                    jobname = f"{field}_{config}_{spw}_{imtype}"

                    match = tbl['JobName'] == jobname
                    if any(match):
                        states = np.unique(tbl['State'][match])
                        if 'RUNNING' in states:
                            jobid = tbl['JobID'][match & (tbl['State'] == 'RUNNING')]
                            log.debug(f"Skipped job {jobname} because it's RUNNING as {set(jobid)}")
                            continue
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

                    # /orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/
                    # group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X130/calibrated/working/
                    # tclean_cube_pars_Sgr_A_st_y_03_TM1_spw27.py

                    basename = f'{field}_{spw}_{imtype}{contsub_suffix}'
                    # basename = "{0}_{1}_spw{2}_{3}".format(field, band, spw, arrayname)

                    tempdir_name = f'{field}_{spw}_{imtype}_{config}_{mousname[6:].replace("/", "_")}'
                    assert any(x in tempdir_name for x in ('7M', 'TM1', 'TP'))

                    # it is safe to remove things beyond here because at this point we're committed
                    # to re-running
                    if '--dry-run' not in sys.argv:
                        if '--remove-failed' in sys.argv:
                            # print(f"Removing files matching '{workdir}/{basename}.*'")
                            failed_files = glob.glob(f'{workdir}/{basename}.*')
                            if any('.image' in x for x in failed_files):
                                print(f"Found a .image in the failed file list: {failed_files}.  Continuing.")
                                # raise ValueError(f"Found a .image in the failed file list: {failed_files}")
                            else:
                                for ff in failed_files:
                                    print(f"Removing {ff}")
                                    shutil.rmtree(ff)

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

                    # DONT start the script from within the appropriate workdir: it puts the log in the wrong place?
                    #os.chdir(f'{workdir}/{tempdir_name}')

                    cmd = (f'/opt/slurm/bin/sbatch --ntasks={ntasks} --cpus-per-task={cpus_per_task} '
                           f'--mem={mem} --output={jobname}_%j.log --job-name={jobname} --account={account} '
                           f'--qos={qos_} --export=ALL --time={jobtime} {runcmd}')

                    if '--dry-run' in sys.argv:
                        if verbose:
                            print(cmd)
                            print()
                        # print(subprocess.check_output('env').decode())
                        # raise
                    else:
                        sbatch = subprocess.check_output(cmd.split())

                        print(f"Started sbatch job with jobid={sbatch.decode()} and parameters {spwpars} and script {scriptname}")

    globals().update(locals())
