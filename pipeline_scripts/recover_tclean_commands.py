"""
Run this script from the base directory 2021.1.00172.L
"""
import os, sys, glob

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

projcode = '2021.1.00172.L'

science_goal_dirs = glob.glob("science_goal*")

def keywriter(x):
    if isinstance(x, str) and 'science_goal' in x:
        return x
    elif isinstance(x, (list, int, float)):
        return x
    else:
        return f'"{x}"'


all_cubepars = {}

for scigoal in science_goal_dirs:
    for group in glob.glob(scigoal+"/*"):
        for member in glob.glob(os.path.join(group, "*")):
            pipeline_run = glob.glob(os.path.join(member, "calibrated/working/pipeline*/html/stage*"))

            cubepars = []
            for stage in pipeline_run:
                logfile = os.path.join(stage, 'casapy.log')
                if os.path.exists(logfile):
                    with open(logfile, 'r') as fh:
                        for line in fh.readlines():
                            if 'tclean( vis' in line:
                                tcleanpars = line.split("tclean")[-1]
                                tcleanpars = eval(f'dict{tcleanpars}')
                                if tcleanpars['specmode'] == 'cube':
                                    cubepars.append(tcleanpars)
                    if cubepars:
                        break

            # only keep the 'iter1' examples
            cubepars = [pars for pars in cubepars
                    if 'iter1' in pars['imagename']]

            all_cubepars[os.path.basename(member)] = cubepars

with open('default_tclean_commands.py', 'w') as tcfh:
    json.dump(all_cubepars, tcfh, indent=2)
