echo "Use job_runner instead"
exit
#!/bin/bash
#SBATCH --job-name=run_reimage_mpi      # Job name
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --ntasks=32                     # 
#SBATCH --mem=256gb                     # Job memory request
#SBATCH --nodes=1                     # Job memory request
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/run_reimage_mpi_%j.log   # Standard output and error log
#SBATCH --qos=adamginsburg-b
#SBATCH --account=adamginsburg

pwd; hostname; date
echo "Memory=${MEM}"

#module load cuda/11.0.207
module load intel/2020.0.166
module load openmpi/4.1.1
#module load libfuse/3.10.4
# test

LOG_DIR=/blue/adamginsburg/adamginsburg/ACES/logs
export LOGFILENAME="${LOG_DIR}/casa_log_mpi_reimaging_${SLURM_JOB_ID}_$(date +%Y-%m-%d_%H_%M_%S).log"

WORK_DIR='/orange/adamginsburg/ACES/rawdata/2021.1.00172.L'
#cd ${WORK_DIR}
# this directory should contain a folder pipeline_scripts/ if any overloaded pipeline scripts are expected
export ACES_ROOTDIR="/orange/adamginsburg/ACES/reduction_ACES/"

# 6.2.1 is the 2021 release version, but it contains the readBlock error
CASAVERSION=casa-6.2.1-7-pipeline-2021.2.0.128
CASAVERSION=casa-6.4.3-2-pipeline-2021.3.0.17
CASAVERSION=casa-6.4.4-31-py3.8
CASAVERSION=casa-6.4.1-12-pipeline-2022.2.0.68
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

cp ~/.casa/config.py ${TMPDIR}

# this appears not to work, but it should
export APPIMAGE_EXTRACT_AND_RUN=1

if [ -z $SLURM_NTASKS ]; then
    echo "FAILURE: SLURM_NTASKS was not specified"
    exit 1
fi

export ACES_DATADIR=/orange/adamginsburg/ACES/rawdata/
export TRYDROPTARGET=True
export SCRIPTLIST='reimaging_script_list.txt' # written in write_tclean_scripts.py
REIMAGING=/orange/adamginsburg/ACES/reduction_ACES/imaging/write_tclean_scripts.py

env

# "dummy" run to print what we'll do
DUMMYRUN=True RUNONCE=False ${MPICASA} -n $SLURM_NTASKS xvfb-run -d ${CASA} --logfile=${LOGFILENAME} --nogui --nologger --rcdir=${TMPDIR} -c "execfile('${REIMAGING}')"

echo "Completed dummy run"

cat $SCRIPTLIST
firstscript=$(head -n 1 $SCRIPTLIST)

if [ -z $firstscript ]; then
    echo "FAILURE: no $firstscript"
    exit 1
fi

echo "Running the following tclean command: $firstscript"
cat $firstscript
echo
echo

cd /blue/adamginsburg/adamginsburg/ACES/workdir 

#echo "Renaming job"
#rename the job to explain what is being cleaned
jn=$(basename $firstscript)
jn=${jn/.py/}
jn=${jn/tclean_cube_pars_Sgr_A_st_/}
#scontrol update jobname=reimaging_${jn}
#echo "Job renamed to reimaging_${jn}"

echo ${MPICASA} -n $SLURM_NTASKS xvfb-run -d ${CASA} --logfile=${LOGFILENAME} --nogui --nologger --rcdir=${TMPDIR} -c "execfile('$firstscript')"
${MPICASA} -n $SLURM_NTASKS xvfb-run -d ${CASA} --logfile=${LOGFILENAME} --nogui --nologger --rcdir=${TMPDIR} -c "execfile('$firstscript')"
