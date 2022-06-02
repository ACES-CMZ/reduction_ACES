#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --nodes=1
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/run_simple_%j.log   # Standard output and error log

env
pwd; hostname; date
echo "Memory=${MEM}"
echo "job name = $jobname"

#module load cuda/11.0.207
module load intel/2020.0.166 openmpi/4.1.1 libfuse/3.10.4

LOG_DIR=/blue/adamginsburg/adamginsburg/ACES/logs
export LOGFILENAME="${LOG_DIR}/casa_log_${jobname}_${SLURM_JOB_ID}_$(date +%Y-%m-%d_%H_%M_%S).log"

WORK_DIR='/blue/adamginsburg/adamginsburg/ACES/workdir/'
cd ${WORK_DIR}
# this directory should contain a folder pipeline_scripts/ if any overloaded pipeline scripts are expected
export ACES_ROOTDIR="/orange/adamginsburg/ACES/reduction_ACES/"

# 6.2.1 is the 2021 release version, but it contains the readBlock error
CASAVERSION=casa-6.4.3-2-pipeline-2021.3.0.17
export CASAPATH=/orange/adamginsburg/casa/${CASAVERSION}
export MPICASA=${CASAPATH}/bin/mpicasa
export CASA=${CASAPATH}/bin/casa
export casapython=${CASAPATH}/lib/py/bin/python3

# CASA MPIrun setup
# we use --export=ALL in srun to export these so CASA knows it's being run in an mpi-enabled environment
export OMPI_COMM_WORLD_SIZE=$SLURM_NTASKS
if [ -z $OMP_NUM_THREADS ]; then
    export OMP_NUM_THREADS=1
fi
export INTERACTIVE=0
export LD_LIBRARY_PATH=${CASAPATH}/lib/:$LD_LIBRARY_PATH
export OPAL_PREFIX="${CASAPATH}/lib/mpi"

export IPYTHONDIR=$SLURM_TMPDIR
export IPYTHON_DIR=$IPYTHONDIR
cp ~/.casa/config.py $SLURM_TMPDIR

if [ -z $SLURM_NTASKS ]; then
    echo "FAILURE: SLURM_NTASKS was not specified"
    exit 1
fi

# try reducing the number of threads by 1 assuming that the 'master' thread can't handle it
# June 1, 2022: try to reduce to ntasks/2 + 1 because this is what pipelines appear to use
mpi_ntasks=$(python -c "print(int(${SLURM_NTASKS} / 2 + 1))")

echo ""
echo "env just before running the command:"
env

${MPICASA} -n ${mpi_ntasks} xvfb-run -d ${CASA} --logfile=${LOGFILENAME} --nogui --nologger --rcdir=$SLURM_TMPDIR -c "execfile('$SCRIPTNAME')"

# the following should not be necessary (it's actually totally redundant) but
# mpicasa has been very severely failing yet reporting COMPLETED instead of
# FAILED
ppid="$!"
wait $ppid
exitcode=$?
exit $exitcode
