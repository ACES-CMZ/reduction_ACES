from parse_weblog import (weblog_names, make_links, get_human_readable_name,
                          get_mous_to_sb_mapping, get_all_fluxes, fluxes_to_table)
import glob
import os
import json

def main():
    if 'weblog_dir' not in locals():
        weblog_dir = os.getenv('WEBLOG_DIR')
        if weblog_dir is None:
            raise ValueError("Set an environmental variable 'WEBLOG_DIR' or a local variable `weblog_dir` to specify where to extract the weblogs.")

    os.chdir(weblog_dir)

    print("Getting MOUS->SOUS mapping")
    mapping = get_mous_to_sb_mapping('2021.1.00172.L')

    if not os.path.exists('humanreadable'):
        os.mkdir('humanreadable')
    weblogs = glob.glob("pipeline*")

    print("Inferring weblog name mappings")
    weblog_maps = weblog_names(weblogs, mapping)

    print("Making links")
    make_links(weblog_maps)

    print("Extracting calibrator fluxes from weblogs")
    fluxes = get_all_fluxes(weblogs)

    print("Dumping fluxes to fluxes.json and fluxes.ipac")
    with open('fluxes.json', 'w') as fh:
        json.dump(fluxes, fh)

    fluxtbl = fluxes_to_table(fluxes)
    for colname in fluxtbl.colnames:
        fluxtbl.rename_column(colname, colname.replace(" ","_"))
    fluxtbl.write('fluxes.ipac', format='ascii.ipac', overwrite=True)
