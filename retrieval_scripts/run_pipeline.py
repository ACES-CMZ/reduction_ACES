"""
Pipeline Running Script
-----------------------

This script is intended to run in the 2021.1.00172.L directory, i.e., the
parent directory of all sous/gous/mous/<project stuff> subdirectories extracted
from the ALMA tarballs.

The script will traverse all subdirectories searching for `scriptForPI.py` files
and will run them if no corresponding `calibrated` directory exists.  If the
calibrated directory exists, the pipeline run will be skipped.

Environmental Variables read by this script:
    ACES_ROOTDIR : the directory containing any scripts that need to be imported.
    This script allows scriptForPIs to be overridden by putting scriptForCalibration
    files in directories following the naming schemem:
        {rootdir}/pipeline_scripts/{sdm}.ms.scriptForCalibration.py
    LOGFILENAME : A CASA log filename to use instead of the default

    RUNONCE: A flag to set - if this is set, it will only run at most
    one scriptForPI and then will quit.  This can be useful if you need to run
    scripts on machines with automatic timeouts, since you don't want a run to
    end partway along.

"""
import os
import sys
import runpy
import glob

try:
    from taskinit import casalog
except ImportError:
    from casatasks import casalog

def logprint(string, origin='almaimf_metadata',
             priority='INFO'):
    print(string)
    casalog.post(string, origin=origin, priority=priority)

if os.getenv('ACES_ROOTDIR') is None:
    try:
        os.environ['ACES_ROOTDIR'] = os.path.split(metadata_tools.__file__)[0]
    except ImportError:
        raise ValueError("metadata_tools not found on path; make sure to "
                         "specify ACES_ROOTDIR environment variable "
                         "or your PYTHONPATH variable to include the directory"
                         " containing the ACES code.")
else:
    sys.path.append(os.getenv('ACES_ROOTDIR'))

rootdir = os.environ['ACES_ROOTDIR']

import shutil
import os

if os.getenv('LOGFILENAME'):
    if os.getenv('LOGFILENAME').startswith('/'):
        casalog.setlogfile(os.getenv('LOGFILENAME'))
    else:
        casalog.setlogfile(os.path.join(os.getcwd(), os.getenv('LOGFILENAME')))
print("CASA log file name is {0}".format(casalog.logfile()))


runonce = bool(os.environ.get('RUNONCE'))

# to do diagnostic plotting, we need the MS, not just the science-only calibrated MS
SPACESAVING = 1
DOSPLIT = True

imaging_script = 'imaging_pipeline_rerun.py'

# only walk the science goal directories
science_goal_dirs = glob.glob("science_goal*")

for scigoal in science_goal_dirs:
    for group in glob.glob(scigoal+"/*"):
        for member in glob.glob(os.path.join(group, "*")):
            dirpath = member
            scriptsforpi = glob.glob(os.path.join(dirpath, "script/*scriptForPI.py"))

            if len(scriptsforpi) == 1:
                scriptforpi = scriptsforpi[0]
            elif len(scriptsforpi) > 1:
                raise ValueError("Too many scripts for PI in {0}".format(dirpath))
            elif len(scriptsforpi) == 0:
                logprint("Skipping directory {0} because it has no scriptForPI"
                         .format(dirpath))
                continue

            curdir = os.getcwd()

            if os.path.exists(os.path.join(dirpath, 'calibrated')):
                if os.path.exists(os.path.join(dirpath, 'calibrated', imaging_script)):
                    logprint("Skipping script {0} in {1} because calibrated "
                             "exists".format(scriptforpi, dirpath), origin='pipeline_runner')
                else:
                    os.chdir(os.path.join(dirpath, 'script'))
                    scriptpath = ("{rootdir}/pipeline_scripts/{imaging_script}"
                                  .format(rootdir=rootdir, imaging_script=imaging_script))

                    shutil.copy(scriptpath, '.')

                    logprint("Running script {0} in {1}".format(imaging_script, dirpath),
                             origin='pipeline_runner')

                    result = runpy.run_path(imaging_script, init_globals=globals())

                    logprint("Done running script {0} in {1}".format(imaging_script, dirpath),
                             origin='pipeline_runner')

                    os.chdir(curdir)
                    
            elif os.path.exists(os.path.join(dirpath, 'calibration')):
                os.chdir(os.path.join(dirpath, 'script'))

                # check for custom scripts
                sdms = glob.glob(os.path.join("../raw/*.asdm.sdm"))
                # reset this each loop so we can search for the custom version
                local_scriptforPI = None
                for sdmfn in sdms:
                    sdm = os.path.split(sdmfn)[-1].split(".")[0]
                    # custom version has to follow this precise name scheme
                    scriptpath = ("{rootdir}/pipeline_scripts/{sdm}.ms.scriptForCalibration.py"
                                  .format(rootdir=rootdir, sdm=sdm))
                    if os.path.exists(scriptpath):
                        shutil.copy(scriptpath, '.')
                        local_scriptforPI = os.path.split(scriptpath)[-1]

                if local_scriptforPI is None:
                    local_scriptforPI = os.path.basename(scriptforpi)

                logprint("Running script {0} in {1}".format(local_scriptforPI, dirpath),
                         origin='pipeline_runner')

                result = runpy.run_path(local_scriptforPI, init_globals=globals())
                # too verbose, includes a ton of junk
                #logprint("result = {0}".format(result),
                #         origin='pipeline_runner')

                logprint("Done running script {0} in {1}".format(local_scriptforPI, dirpath),
                         origin='pipeline_runner')

                os.chdir(curdir)

                if runonce:
                    sys.exit(0)
            else:
                raise ValueError("Landed in the wrong directory.")

logprint("Completed run_pipeline.py")
