#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --ntasks=6
#SBATCH --mem=24gb # Job memory request PER NODE
#SBATCH --nodes=1 # exactly 1 node
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --qos=adamginsburg-b
#SBATCH --account=adamginsburg
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/ACES_retrieval_%j.log
#SBATCH --job-name=ACES_retrieval
#SBATCH --export=ALL


date

. ~/.gh_token
echo $GITHUB_TOKEN

cd /orange/adamginsburg/ACES/rawdata

echo "test import"
/orange/adamginsburg/miniconda3/envs/python39/bin/python -c "import zipfile"
echo "Retrieve data"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/retrieval_scripts/retrieve_data.py keflavich True
echo "Retrieve weblogs"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/retrieval_scripts/retrieve_weblogs.py keflavich


export WEBLOG_DIR=/orange/adamginsburg/web/secure/ACES/weblogs/
export WEBLOG_DIR='/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/weblogs'

echo "Make links"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/retrieval_scripts/make_humanreadable_links.py
echo "Update github"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/hipergator_scripts/ghapi_update.py


echo "Make 7m mosaic"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/imaging/mosaic_7m.py
echo "Make 12m mosaic"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/imaging/mosaic_12m.py
echo "Make TP mosaic"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/imaging/mosaic_TP.py

# technically shouldn't need to be re-run, but as I add new mosaics, it will
# ln -s /orange/adamginsburg/ACES/mosaics/*png /orange/adamginsburg/web/secure/ACES/mosaics/
