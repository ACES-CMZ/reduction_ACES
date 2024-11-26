export USE_DASK=False
export USE_LOCAL=True

for MOLNAME in CH3CHO HC15N SO21 H13CN HN13C H13COp CS21 HC3N HCOP SiO21 SO32 HNCO_7m12mTP HNCO NSplus H40a; do

    if [ -e /orange/adamginsburg/ACES/mosaics/cubes/${MOLNAME}_CubeMosaic.fits ]; then

        export MOLNAME

        echo "Giant ${MOLNAME} cube"
        jobid=$(DO_PV=True DOWNSAMPLE=True sbatch --job-name=${MOLNAME}_giantcube_analysis --output=/red/adamginsburg/ACES/logs/ACES_${MOLNAME}_analysis_%j.log --export=ALL --account=astronomy-dept --qos=astronomy-dept-b --ntasks=32 --nodes=1 --mem=256gb --time=96:00:00 --wrap "/orange/adamginsburg/miniconda3/envs/python312/bin/python /orange/adamginsburg/ACES/reduction_ACES/aces/analysis/giantcube_downsample.py")
    else
        echo "${MOLNAME}_CubeMosaic.fits does not exist"
        ls -lh /orange/adamginsburg/ACES/mosaics/cubes/${MOLNAME}*fits
    fi
done
