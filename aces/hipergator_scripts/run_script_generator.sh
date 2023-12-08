#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=4gb                     # Job memory request
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/make_ACES_scripts_%j.log
#SBATCH --export=ALL
#SBATCH --job-name=make_ACES_scripts
#SBATCH --qos=adamginsburg-b
#SBATCH --account=adamginsburg
pwd; hostname; date

export ACES="/orange/adamginsburg/ACES/"
export WORK_DIR="/orange/adamginsburg/ACES/reduction_ACES/"
export WORK_DIR="/blue/adamginsburg/adamginsburg/ACES/workdir"
export ACES_DATADIR="${ACES}/data/"
export WEBLOG_DIR="${ACES}/data/2021.1.00172.L/weblogs/"

module load git

which python
which git

git --version
echo $?

export PYPATH=/orange/adamginsburg/miniconda3/envs/python39/bin/
export IPYTHON=/orange/adamginsburg/miniconda3/envs/python39/bin/ipython

cd ${WORK_DIR}
echo ${WORK_DIR}

export ACES_ROOTDIR="/orange/adamginsburg/ACES/reduction_ACES/"
export SCRIPT_DIR="${ACES_ROOTDIR}/analysis"
export PYTHONPATH=$SCRIPT_DIR

echo $LOGFILENAME

export NO_PROGRESSBAR='True'
export ENVIRON='BATCH'

# TODO: get this from ntasks?
#export DASK_THREADS=8
export DASK_THREADS=$SLURM_NTASKS

export TRYDROPTARGET=True
export TEMPORARY_WORKING_DIRECTORY="/blue/adamginsburg/adamginsburg/ACES/workdir"

env

echo "Recovering tclean commands"
${PYPATH}/aces_recover_tclean_commands || exit 1
echo "Writing tclean scripts"
${PYPATH}/aces_write_tclean_scripts || exit 1
echo "Calling job running"
${PYPATH}/aces_job_runner --verbose=True --parallel || exit 1
echo "Linking repipeline weblogs"
${PYPATH}/aces_link_repipeline_weblogs || exit 1
echo "Done!"

# || exit 1 should quit if any task fails
#${IPYTHON} ${ACES_ROOTDIR}/pipeline_scripts/recover_tclean_commands.py || exit 1
#${IPYTHON} ${ACES_ROOTDIR}/imaging/write_tclean_scripts.py || exit 1
#${IPYTHON} ${ACES_ROOTDIR}/hipergator_scripts/job_runner.py --verbose=True || exit 1
#${IPYTHON} ${ACES_ROOTDIR}/hipergator_scripts/link_repipeline_weblogs.py || exit 1
