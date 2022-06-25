"""
Run this with some care; it will require manual pushing
"""

import json

with open(f"{rootdir}/pipeline_scripts/default_tclean_commands.json", "r") as fh:
    default_commands = json.load(fh)

with open(f"{rootdir}/pipeline_scripts/override_tclean_commands.json", "r") as fh:
    override_commands = json.load(fh)

from merge_tclean_commands import commands

for key in commands:
    if 'TM' in key:
        if 'spw33' not in commands[key]['tclean_cube_pars']:
            spw33pars = commands[key]['tclean_cube_pars']['spw35']
            spw33pars['nchan'] = 3836
            spw33pars['start'] = "97.6660537907GHz"
            spw33pars['width'] = "0.4882381MHz"
            spw33pars['threshold'] = '0.01Jy'
            override_commands[key]['tclean_cube_pars'] = spw33pars
        elif commands[key]['tclean_cube_pars']['spw33']['nchan'] < 3800:
            spw33pars = {}
            spw33pars['nchan'] = 3836
            spw33pars['start'] = "97.6660537907GHz"
            spw33pars['width'] = "0.4882381MHz"
            spw33pars['threshold'] = '0.01Jy'
            override_commands[key]['tclean_cube_pars'] = spw33pars

with open(f"{rootdir}/pipeline_scripts/override_tclean_commands.json", "w") as fh:
    json.dump(fh, override_commands)
