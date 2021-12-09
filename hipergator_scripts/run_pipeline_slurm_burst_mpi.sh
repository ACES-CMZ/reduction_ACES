#!/bin/bash
#SBATCH --job-name=run_pipeline_mpi      # Job name
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --ntasks=8                     # try w/8 to improve queue time?
#SBATCH --mem=32gb                     # Job memory request
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --output=run_pipeline_mpi_%j.log   # Standard output and error log
#SBATCH --qos=adamginsburg-b
#SBATCH --account=adamginsburg


module load intel/2019.1.144 openmpi/4.0.1

WORK_DIR='/orange/adamginsburg/ACES/rawdata/2021.1.00172.L'
cd ${WORK_DIR}
# this directory should contain a folder pipeline_scripts/ if any overloaded pipeline scripts are expected
export ACES_ROOTDIR="/orange/adamginsburg/ACES/reduction_ACES/"

CASAVERSION=casa-6.2.1-7-pipeline-2021.2.0.128
export MPICASA=/orange/adamginsburg/casa/${CASAVERSION}/bin/mpicasa
export CASA=/orange/adamginsburg/casa/${CASAVERSION}/bin/casa

# since we're bursting, be careful not to partially start a pipeline run
export RUNONCE=True

# casa's python requires a DISPLAY for matplot so create a virtual X server
xvfb-run -d ${MPICASA} -n 8 ${CASA} --pipeline --nogui --nologger -c "execfile('$ACES_ROOTDIR/run_pipeline.py')"
