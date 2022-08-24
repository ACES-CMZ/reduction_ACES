"""
This will create (or overwrite) existing tclean commands with label 'aggregate_high' and therefore should not be run automatically.
"""

import json
import os
from astropy import log

from aces.pipeline_scripts.merge_tclean_commands import commands
from aces import conf
rootdir = os.path.join(conf.basepath, "reduction_ACES")

if __name__ == "__main__":

    if os.getenv('ACES_ROOTDIR') is not None:
        log.warning(f"Overridding default rootdir={rootdir} with rootdir={os.environ['ACES_ROOTDIR']}")
        rootdir = os.environ['ACES_ROOTDIR']

    with open(f"{rootdir}/aces/pipeline_scripts/default_tclean_commands.json", "r") as fh:
        default_commands = json.load(fh)

    override_file = f"{rootdir}/aces/pipeline_scripts/override_tclean_commands.json"
    with open(override_file, "r") as fh:
        override_commands = json.load(fh)

    ncmds = (len(override_commands))

    log.debug(f"There are {len(override_commands)} override commands and {len(commands)} commands.")

    for key in commands:
        if 'TM' in key:
            if key not in override_commands:
                override_commands[key] = {}
            if 'tclean_cont_pars' not in override_commands[key]:
                override_commands[key]['tclean_cont_pars'] = {}
            if 'aggregate_high' not in commands[key]['tclean_cont_pars']:
                print(f"Adding {key}")
                pars = commands[key]['tclean_cont_pars']['aggregate']
                pars['imagename'] = pars['imagename'].replace('25_27_29_31_33_35', '33_35')

                # split all the spw selections such that only 33 and 35 are kept
                spwsel = pars['spw']
                spwsel = ["33:" + x.split("33:")[-1] for x in spwsel]
                pars['spw'] = spwsel

                override_commands[key]['tclean_cont_pars']['aggregate_high'] = pars

    assert len(override_commands) >= ncmds

    log.debug(f"Overwriting to {override_file}")
    with open(override_file, "w") as fh:
        json.dump(override_commands, fh, indent=2)