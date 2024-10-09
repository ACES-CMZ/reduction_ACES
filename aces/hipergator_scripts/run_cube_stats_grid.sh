#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --ntasks=64
#SBATCH --nodes=1
#SBATCH --mem=256gb                     # Job memory request
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --output=/red/adamginsburg/ACES/logs/cube_stats_grid_ACES_%j.log
#SBATCH --export=ALL
#SBATCH --job-name=cube_stats_grid_ACES
#SBATCH --qos=astronomy-dept-b
#SBATCH --account=astronomy-dept
pwd; hostname; date

export WORK_DIR="/blue/adamginsburg/adamginsburg/ACES/workdir"

module load git

which python
which git

git --version
echo $?

export IPYTHON=/orange/adamginsburg/miniconda3/envs/python310/bin/ipython

cd ${WORK_DIR}
echo ${WORK_DIR}

export ACES_ROOTDIR="/orange/adamginsburg/ACES/reduction_ACES/"
export SCRIPT_DIR="${ACES_ROOTDIR}/analysis"
export PYTHONPATH=$SCRIPT_DIR

echo $LOGFILENAME

export NO_PROGRESSBAR='True'
export ENVIRON='BATCH'
export JOBNAME=cube_stats_grid_ACES
export jobname=$JOBNAME

# TODO: get this from ntasks?
#export DASK_THREADS=8
export DASK_THREADS=$SLURM_NTASKS

#export START_FROM_CACHED=False

env

echo "Starting cube stats grid feathered"
/orange/adamginsburg/miniconda3/envs/python310/bin/aces_cube_stats_grid_feathered
echo "Done with cube stats grid feathered"

/orange/adamginsburg/miniconda3/envs/python310/bin/aces_cube_stats_grid
