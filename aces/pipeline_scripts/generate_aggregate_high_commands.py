"""
This will create (or overwrite) existing tclean commands with label 'aggregate_high' and therefore should not be run automatically.

(actually this should run automatically)
"""

import json
import os
from astropy import log

from aces import conf
rootdir = os.path.join(conf.basepath, "reduction_ACES")

if os.getenv('ACES_ROOTDIR') is not None:
    log.warning(f"Overridding default rootdir={rootdir} with rootdir={os.environ['ACES_ROOTDIR']}")
    rootdir = os.environ['ACES_ROOTDIR']

pipedir = os.path.dirname(__file__)


def main():

    with open(f"{pipedir}/default_tclean_commands.json", "r") as fh:
        commands = json.load(fh)

    aggregate_high_file = f"{pipedir}/aggregate_high_tclean_commands.json"
    if os.path.exists(aggregate_high_file):
        with open(aggregate_high_file, "r") as fh:
            aggregate_high_commands = json.load(fh)
    else:
        aggregate_high_commands = {}

    ncmds = (len(aggregate_high_commands))

    log.debug(f"There are {len(aggregate_high_commands)} aggregate_high commands and {len(commands)} commands.")

    for key in commands:
        if 'TM' in key:
            if key not in aggregate_high_commands:
                aggregate_high_commands[key] = {}
            if 'tclean_cont_pars' not in aggregate_high_commands[key]:
                aggregate_high_commands[key]['tclean_cont_pars'] = {}
            if 'aggregate_high' not in commands[key]['tclean_cont_pars']:
                print(f"Adding {key}")
                pars = commands[key]['tclean_cont_pars']['aggregate']
                pars['imagename'] = pars['imagename'].replace('25_27_29_31_33_35', '33_35')

                # split all the spw selections such that only 33 and 35 are kept
                spwsel = pars['spw']
                spwsel = ["33:" + x.split("33:")[-1] for x in spwsel]
                pars['spw'] = spwsel

                aggregate_high_commands[key]['tclean_cont_pars']['aggregate_high'] = pars

    assert len(aggregate_high_commands) >= ncmds

    log.debug(f"Overwriting to {aggregate_high_file}")
    with open(aggregate_high_file, "w") as fh:
        json.dump(aggregate_high_commands, fh, indent=2)
