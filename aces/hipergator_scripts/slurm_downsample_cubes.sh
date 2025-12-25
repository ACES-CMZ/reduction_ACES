export USE_DASK=True
# used by giantcube_cuts but not by giantcube_downsample export USE_LOCAL=True

for MOLNAME in CH3CHO HC15N SO21 H13CN HN13C H13COp CS21 HC3N HCOP SiO21 SO32 HNCO_7m12mTP HNCO NSplus H40a; do

    if [ -e /orange/adamginsburg/ACES/mosaics/cubes/${MOLNAME}_CubeMosaic.fits ]; then

        export MOLNAME

        echo "Giant ${MOLNAME} cube"
        jobid=$(DO_PV=True DOWNSAMPLE=True sbatch \
            --job-name=${MOLNAME}_giantcube_downsample --output=/blue/adamginsburg/adamginsburg/ACES/logs/ACES_${MOLNAME}_downsample_%j.log --export=ALL --account=astronomy-dept \
            --qos=astronomy-dept-b --ntasks=16 --nodes=1 --mem=128gb --time=96:00:00 \
            --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python312/bin/python /orange/adamginsburg/ACES/reduction_ACES/aces/analysis/giantcube_downsample.py")
    else
        echo "${MOLNAME}_CubeMosaic.fits does not exist"
        ls -lh /orange/adamginsburg/ACES/mosaics/cubes/${MOLNAME}*fits
    fi
done
