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
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/ACES_CS21_analysis_%j.log
#SBATCH --job-name=ACES_CS21_analysis
#SBATCH --export=ALL

date

. ~/.gh_token
echo $GITHUB_TOKEN

cd /blue/adamginsburg/adamginsburg/ACES/workdir/
pwd

export USE_DASK=True

if [ -e /orange/adamginsburg/ACES/mosaics/cubes/CS21_CubeMosaic.fits ]; then

    echo "test import"
    /orange/adamginsburg/miniconda3/envs/python39/bin/python -c "import zipfile" || exit 1

    echo "Giant CS 21 cube"
    /orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/aces/analysis/giantcube_cuts.py || exit 1
else
    echo "CS21_CubeMosaic.fits does not exist"
    ls -lh CS*fits
fi
