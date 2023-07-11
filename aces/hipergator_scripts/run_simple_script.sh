#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --nodes=1
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/run_simple_%j.log   # Standard output and error log

# https://stackoverflow.com/a/2853811/814354
set -o xtrace

env
pwd; hostname; date
echo "Memory=${MEM}"
echo "job name = $jobname"

if [ -z $jobname ]; then
    export jobname=$SLURM_JOB_NAME
    echo "job name = $jobname (set from SLURM_JOB_NAME=${SLURM_JOB_NAME}"
fi

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
# CASA 6.4.3-2 does not seem to be a real version anywhere, so we move now to 6.4.1-12
CASAVERSION=casa-6.4.1-12-pipeline-2022.2.0.68
# 6.4.1-12 above does not work - it just gives OOM errors for any operation
CASAVERSION=casa-6.5.5-21-py3.8
CASAVERSion=casa-6.5.7-1-py3.8.el7
# none of the new versions work, maybe this one will?
# Does not support nmajor? CASAVERSION=casa-6.2.1-7-pipeline-2021.2.0.128
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
echo "mpi_ntasks=${mpi_ntasks}"
pwd

if [ ${mpi_ntasks} -gt 1 ]; then
    echo ${MPICASA} -n ${mpi_ntasks} xvfb-run -d ${CASA} --logfile=${LOGFILENAME} --nogui --nologger --rcdir=$SLURM_TMPDIR -c "execfile('$SCRIPTNAME')"
    ${MPICASA} -n ${mpi_ntasks} xvfb-run -d ${CASA} --logfile=${LOGFILENAME} --nogui --nologger --rcdir=$SLURM_TMPDIR -c "execfile('$SCRIPTNAME')" || exit 1
else
    echo xvfb-run -d ${CASA} --logfile=${LOGFILENAME} --nogui --nologger --rcdir=$SLURM_TMPDIR -c "execfile('$SCRIPTNAME')"
    xvfb-run -d ${CASA} --logfile=${LOGFILENAME} --nogui --nologger --rcdir=$SLURM_TMPDIR -c "execfile('$SCRIPTNAME')" || exit 1
fi

# the following should not be necessary (it's actually totally redundant) but
# mpicasa has been very severely failing yet reporting COMPLETED instead of
# FAILED
# ppid="$!"
# wait $ppid
# exitcode=$?
# echo "pid=$ppid exitcode=$exitcode"
# exit $exitcode
