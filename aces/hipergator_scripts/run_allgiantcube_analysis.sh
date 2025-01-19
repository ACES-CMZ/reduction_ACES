#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --ntasks=64
#SBATCH --mem=256gb
#SBATCH --nodes=1 # exactly 1 node
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --qos=astronomy-dept-b
#SBATCH --account=astronomy-dept
#SBATCH --output=/red/adamginsburg/ACES/logs/ACES_AllCube_analysis_%j.log
#SBATCH --job-name=ACES_AllCube_analysis
#SBATCH --export=ALL

date

. ~/.gh_token
echo $GITHUB_TOKEN

cd /red/adamginsburg/ACES/workdir/
pwd

# use_dask=False resulted in an immediate OOM, which doesn't make sense
# but at least as of 2025/01/19, non-dask works and dask fails
export USE_DASK=False
export USE_LOCAL=True

# 2025/01/16: attempt a different dask approach
export DASK_CLIENT=threads

#for MOLNAME in HC3N; do
for MOLNAME in CH3CHO NSplus H40a HC15N SO21 H13CN HN13C H13COp CS21 HC3N HCOP SiO21 SO32 HNCO_7m12mTP; do

    if [ -e /orange/adamginsburg/ACES/mosaics/cubes/${MOLNAME}_CubeMosaic.fits ]; then

        #echo "test import"
        #/orange/adamginsburg/miniconda3/envs/python312/bin/python -c "import zipfile" || exit 1

        export MOLNAME
        # optional
        export DOWNSAMPLE=False
        export DO_PV=False

        echo "Giant ${MOLNAME} cube"
        #/red/adamginsburg/miniconda3/envs/python312/bin/python /orange/adamginsburg/ACES/reduction_ACES/aces/analysis/giantcube_cuts.py || exit 1
        sbatch --job-name=aces_giantcube_analysis_${MOLNAME}_dask${DASK_CLIENT} \
               --output=/red/adamginsburg/ACES/logs/aces_giantcube_analysis_${MOLNAME}_dask${DASK_CLIENT}_%j.log \
               --account=astronomy-dept --qos=astronomy-dept-b --ntasks=64 --nodes=1 \
               --mem=256gb --time=96:00:00 \
               --export=MOLNAME=${MOLNAME},USE_DASK=${USE_DASK},USE_LOCAL=${USE_LOCAL},DASK_CLIENT=${DASK_CLIENT},DOWNSAMPLE=${DOWNSAMPLE},DO_PV=${DO_PV} \
               --wrap /red/adamginsburg/miniconda3/envs/python312/bin/aces_giantcube_analysis
    else
        echo "${MOLNAME}_CubeMosaic.fits does not exist"
        ls -lh ${MOLNAME}*fits
    fi
done
