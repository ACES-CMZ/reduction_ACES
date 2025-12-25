date

. ~/.gh_token
echo $GITHUB_TOKEN

cd /blue/adamginsburg/adamginsburg/ACES/workdir/
pwd

# use_dask=False resulted in an immediate OOM, which doesn't make sense
export USE_DASK=True
export USE_LOCAL=True

# 2025/01/16: attempt a different dask approach
export DASK_CLIENT=slurm

for MOLNAME in HC3N; do
#for MOLNAME in CH3CHO NSplus H40a HC15N SO21 H13CN HN13C H13COp CS21 HC3N HCOP SiO21 SO32 HNCO_7m12mTP; do

    if [ -e /orange/adamginsburg/ACES/mosaics/cubes/${MOLNAME}_CubeMosaic.fits ]; then

        #echo "test import"
        #/orange/adamginsburg/miniconda3/envs/python312/bin/python -c "import zipfile" || exit 1

        export MOLNAME
        # optional
        export DOWNSAMPLE=False
        export DO_PV=False

        echo "Giant ${MOLNAME} cube"
        #/blue/adamginsburg/adamginsburg/miniconda3/envs/python312/bin/python /orange/adamginsburg/ACES/reduction_ACES/aces/analysis/giantcube_cuts.py || exit 1
        sbatch --job-name=aces_giantcube_analysis_${MOLNAME}_daskslurm \
               --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_giantcube_analysis_${MOLNAME}_daskslurm_%j.log \
               --account=astronomy-dept --qos=astronomy-dept-b --ntasks=16 --nodes=1 \
               --mem=64gb --time=96:00:00 \
               --export=MOLNAME=${MOLNAME},USE_DASK=${USE_DASK},USE_LOCAL=${USE_LOCAL},DASK_CLIENT=${DASK_CLIENT},DOWNSAMPLE=${DOWNSAMPLE},DO_PV=${DO_PV} \
               --wrap /blue/adamginsburg/adamginsburg/miniconda3/envs/python312/bin/aces_giantcube_analysis
    else
        echo "${MOLNAME}_CubeMosaic.fits does not exist"
        ls -lh ${MOLNAME}*fits
    fi
done

