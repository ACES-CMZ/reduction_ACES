from parse_weblog import (weblog_names, make_links, get_human_readable_name,
                          get_mous_to_sb_mapping, get_all_fluxes, fluxes_to_table)
import glob
import os
import json

mapping = get_mous_to_sb_mapping('2021.1.00172.L')

os.chdir('/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/weblogs')
if not os.path.exists('humanreadable'):
    os.mkdir('humanreadable')
weblogs = glob.glob("pipeline*")

weblog_maps = weblog_names(weblogs, mapping)

make_links(weblog_maps)

fluxes = get_all_fluxes(weblogs)

with open('fluxes.json', 'w') as fh:
    json.dump(fluxes, fh)

fluxtbl = fluxes_to_table(fluxes)
for colname in fluxtbl.colnames:
    fluxtbl.rename_column(colname, colname.replace(" ","_"))
fluxtbl.write('fluxes.ipac', format='ascii.ipac')
