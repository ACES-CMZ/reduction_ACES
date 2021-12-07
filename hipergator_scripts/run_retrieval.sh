#!/bin/bash
#SBATCH --mail-type=NONE          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=adamginsburg@ufl.edu     # Where to send mail
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --ntasks=1
#SBATCH --mem=4gb # Job memory request PER NODE
#SBATCH --nodes=1 # exactly 1 node
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --qos=adamginsburg-b
#SBATCH --account=adamginsburg
#SBATCH --output=ACES_retrieval_%j.log
#SBATCH --job-name=ACES_retrieval
#SBATCH --export=ALL


date

cd /orange/adamginsburg/ACES/rawdata

echo "Retrieve data"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/retrieval_scripts/retrieve_data.py keflavich
echo "Retrieve weblogs"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/retrieval_scripts/retrieve_weblogs.py keflavich


export WEBLOG_DIR=/orange/adamginsburg/web/secure/ACES/weblogs/

echo "Make links"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/retrieval_scripts/make_humanreadable_links.py
echo "Update github"
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/hipergator_scripts/ghapi_update.py
