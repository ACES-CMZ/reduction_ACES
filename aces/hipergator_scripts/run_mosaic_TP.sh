#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --ntasks=32
#SBATCH --mem=256gb
#SBATCH --nodes=1 # exactly 1 node
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --qos=astronomy-dept-b
#SBATCH --account=astronomy-dept
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/ACES_mosaicTP_%j.log
#SBATCH --job-name=ACES_mosaicTP
#SBATCH --export=ALL


date

. ~/.gh_token
echo $GITHUB_TOKEN

cd /orange/adamginsburg/ACES/rawdata

echo "test import"
/orange/adamginsburg/miniconda3/envs/python39/bin/python -c "import zipfile" || exit 1
#echo "Retrieve data"
#/orange/adamginsburg/miniconda3/envs/python39/bin/aces_retrieve_data keflavich True True || exit 1
#echo "Retrieve weblogs"
#/orange/adamginsburg/miniconda3/envs/python39/bin/aces_retrieve_weblogs keflavich || exit 1


#export WEBLOG_DIR=/orange/adamginsburg/web/secure/ACES/weblogs/
#export WEBLOG_DIR='/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/weblogs'
#
#echo "Make links"
#/orange/adamginsburg/miniconda3/envs/python39/bin/aces_make_humanreadable_links || exit 1
#echo "Update github"
#/orange/adamginsburg/miniconda3/envs/python39/bin/aces_ghapi_update || exit 1


#echo "Make 7m mosaic"
#/orange/adamginsburg/miniconda3/envs/python39/bin/aces_mosaic_7m || exit 1
echo "Make TP mosaic"
/orange/adamginsburg/miniconda3/envs/python39/bin/aces_mosaic_TP || exit 1
#echo "Make 12m mosaic"
#/orange/adamginsburg/miniconda3/envs/python39/bin/aces_mosaic_12m || exit 1

# technically shouldn't need to be re-run, but as I add new mosaics, it will
# (but this causes a 'failed' error)
#ln -s /orange/adamginsburg/ACES/mosaics/*png /orange/adamginsburg/web/secure/ACES/mosaics/
