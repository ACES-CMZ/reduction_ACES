#!/bin/bash
#SBATCH --job-name=run_pipeline_mpi      # Job name
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --ntasks=32                     # 
#SBATCH --mem=128gb                     # Job memory request
#SBATCH --nodes=1
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/run_pipeline_mpi_%j.log   # Standard output and error log
#SBATCH --qos=adamginsburg-b
#SBATCH --account=adamginsburg

env
pwd; hostname; date
echo "Memory=${MEM}"

#module load cuda/11.0.207
module load intel/2020.0.166
module load openmpi/4.1.1 
#module load libfuse/3.10.4

LOG_DIR=/blue/adamginsburg/adamginsburg/ACES/logs
export LOGFILENAME="${LOG_DIR}/casa_log_mpi_pipeline_${SLURM_JOB_ID}_$(date +%Y-%m-%d_%H_%M_%S).log"

WORK_DIR='/orange/adamginsburg/ACES/rawdata/2021.1.00172.L'
cd ${WORK_DIR}
# this directory should contain a folder pipeline_scripts/ if any overloaded pipeline scripts are expected
export ACES_ROOTDIR="/orange/adamginsburg/ACES/reduction_ACES/"

CASAVERSION=casa-6.2.1-7-pipeline-2021.2.0.128
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



if [ -z $SLURM_NTASKS ]; then
    echo "FAILURE: SLURM_NTASKS was not specified"
    exit 1
fi

# since we're using a burst node, be careful not to partially start a pipeline run for every folder
export RUNONCE=True

#echo srun a.out
#srun --export=ALL --mpi=pmix_v3 /orange/adamginsburg/test_mpi/a.out
#
#echo srun --export=ALL --mpi=pmix_v3 xvfb-run -d ${CASA} --nogui --nologger -c 'print("TEST")'
#srun --export=ALL --mpi=pmix_v3 xvfb-run -d ${CASA} --nogui --nologger --ipython-dir=/tmp -c 'print("TEST")'
#echo "Done"
#
#
#srun --export=ALL --mpi=pmix_v3 xvfb-run -d ${CASA} --nogui --nologger --ipython-dir=/tmp\
#    -c 'from casampi.MPIEnvironment import MPIEnvironment; from casampi.MPICommandClient import MPICommandClient; print(f"MPI enabled: {MPIEnvironment.is_mpi_enabled}")'
#srun --export=ALL --mpi=pmix_v3 xvfb-run -d ${CASA} --nogui --nologger --ipython-dir=/tmp -c 'import os, pprint; print(pprint.pprint(dict(os.environ)))'


RUNSCRIPTS=False /orange/adamginsburg/casa/${CASAVERSION}/bin/python3 /orange/adamginsburg/ACES/reduction_ACES/retrieval_scripts/run_pipeline.py > /tmp/scriptlist

echo "Script List"
cat /tmp/scriptlist
echo "***********"

for script in $(cat /tmp/scriptlist); do 
    IFS='/' read -r -a array <<< "$script"
    mous=${array[2]}
    echo "MOUS is ${mous}, script is ${script}"
    export LOGFILENAME="${LOG_DIR}/casa_log_mpi_pipeline_${mous}_${SLURM_JOB_ID}_$(date +%Y-%m-%d_%H_%M_%S).log"
    cwd=$(pwd)

    cd $(dirname $script)
    pwd
    echo $script
    echo srun --export=ALL --mpi=pmix_v3 ${CASA} --logfile=${LOGFILENAME} --pipeline --nogui --nologger --ipython-dir=/tmp -c "execfile('${script}')"
    #srun --export=ALL --mpi=pmix_v3 
    echo ${MPICASA} -n $SLURM_NTASKS ${CASA} --logfile=${LOGFILENAME} --pipeline --nogui --nologger --ipython-dir=/tmp -c "execfile('$(basename ${script})')"
    ${MPICASA} -n $SLURM_NTASKS ${CASA} --logfile=${LOGFILENAME} --pipeline --nogui --nologger --ipython-dir=/tmp -c "execfile('$(basename ${script})')"
    cd $cwd
done
echo "Done"
