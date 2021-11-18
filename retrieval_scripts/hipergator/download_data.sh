#!bash

cd /orange/adamginsburg/ACES/rawdata

/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/retrieval_scripts/retrieve_data.py keflavich
/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/retrieval_scripts/retrieve_weblogs.py keflavich

export WEBLOG_DIR='/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/weblogs'

/orange/adamginsburg/miniconda3/envs/python39/bin/python /orange/adamginsburg/ACES/reduction_ACES/retrieval_scripts/make_humanreadable_links.py
