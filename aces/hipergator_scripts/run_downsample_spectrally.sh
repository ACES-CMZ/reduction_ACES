date

. ~/.gh_token
echo $GITHUB_TOKEN

cd /blue/adamginsburg/adamginsburg/ACES/workdir/
pwd


#for MOLNAME in HC15N SO21 H13CN HN13C H40a CH3CHO NSplus; do # H13COp CS21 HC3N HCOP SiO21 HNCO_7m12mTP; do
for MOLNAME in HCOP_mopra HNCO_7m12mTP HCOP_noTP CH3CHO NSplus H40a HC15N SO21 H13CN HN13C H13COp CS21 HC3N HCOP SiO21 SO32; do
#for MOLNAME in HNCO_7m12mTP HCOP_noTP H13CN HCOP; do
#for MOLNAME in HNCO_7m12mTP; do

    if [ -e /orange/adamginsburg/ACES/mosaics/cubes/${MOLNAME}_CubeMosaic.fits ]; then

        #echo "test import"
        #/orange/adamginsburg/miniconda3/envs/python312/bin/python -c "import zipfile" || exit 1

        export MOLNAME
        export USE_DASK=True
        export NUM_CORES=32 # could also rely on slurm for this

        echo "Giant ${MOLNAME} cube - spectrally downsampling"
        #/blue/adamginsburg/adamginsburg/miniconda3/envs/python312/bin/python /orange/adamginsburg/ACES/reduction_ACES/aces/analysis/giantcube_cuts.py || exit 1
        sbatch --job-name=aces_downsample_spectrally_${MOLNAME} \
               --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_downsample_spectrally_${MOLNAME}_%j.log \
               --account=astronomy-dept --qos=astronomy-dept-b --ntasks=${NUM_CORES} --nodes=1 \
               --cpus-per-task=1 \
               --mem=256gb --time=96:00:00 \
               --export=MOLNAME=${MOLNAME},USE_DASK=${USE_DASK},NUM_CORES=${NUM_CORES} \
               --wrap '/blue/adamginsburg/adamginsburg/miniconda3/envs/python312/bin/python -c "from aces.analysis import downsample_all_spectrally; downsample_all_spectrally.main()"'
    else
        echo "${MOLNAME}_CubeMosaic.fits does not exist"
        ls -lh ${MOLNAME}*fits
    fi
done

# one=line version
# export MOLNAME=SiO21 DOWNSAMPLE=True DOPV=True mem=256 USE_LOCAL=True USE_DASK=False;            sbatch --job-name=aces_giantcube_analysis_${MOLNAME}${DASK}${DASK_CLIENT}                --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_giantcube_analysis_${MOLNAME}${DASK}${DASK_CLIENT}_%j.log                --account=astronomy-dept --qos=astronomy-dept-b --ntasks=32 --nodes=1                --mem=${mem}gb --time=96:00:00                --export=MOLNAME=${MOLNAME},USE_DASK=${USE_DASK},USE_LOCAL=${USE_LOCAL},DASK_CLIENT=${DASK_CLIENT},DOWNSAMPLE=${DOWNSAMPLE},DO_PV=${DO_PV}                --wrap /red/adamginsburg/miniconda3/envs/python312/bin/aces_giantcube_analysis
