import json
import itertools
from astropy import log
from astropy.table import Table
import os
from aces import conf

# old version
# if os.getenv('ACES_ROOTDIR') is None:
#     raise ValueError("Specify ACES_ROOTDIR environment variable ")
# else:
#     rootdir = os.environ['ACES_ROOTDIR']

if 'verbose' not in locals():
    verbose = False

pipedir = os.path.dirname(__file__)


def merge_aggregate(commands):
    from aces.pipeline_scripts import generate_aggregate_high_commands
    generate_aggregate_high_commands.main()

    with open(f"{pipedir}/aggregate_high_tclean_commands.json", "r") as fh:
        aggregate_high_commands = json.load(fh)

    for sbname, allpars in aggregate_high_commands.items():
        if sbname not in commands:
            log.warning(f"SB {sbname} was not in the default tclean commands; using aggregate before override")
            commands[sbname] = allpars
            continue
        for partype, replacements in allpars.items():
            if partype in ('manual QA', ):
                continue
            for spwsel, tcpars in replacements.items():
                if spwsel in commands[sbname][partype]:
                    # we're replacing/overriding arguments here
                    for key, val in tcpars.items():
                        orig = commands[sbname][partype][spwsel][key]
                        if verbose:
                            print(f"{sbname} {partype} {spwsel}: {key}: {orig} -> {val}")
                        commands[sbname][partype][spwsel][key] = val
                else:
                    # but if a spw was totally skipped, we replace it
                    # with the override version
                    if verbose:
                        print(f"SPW {spwsel} was skipped in {sbname} and is being replaced")
                    commands[sbname][partype][spwsel] = tcpars
    return commands


def merge_override(commands):
    with open(f"{pipedir}/override_tclean_commands.json", "r") as fh:
        override_commands = json.load(fh)

    for sbname, allpars in override_commands.items():
        if sbname not in commands:
            log.warning(f"SB {sbname} was not in the default tclean commands; using override only")
            commands[sbname] = allpars
            continue
        for partype, replacements in allpars.items():
            if partype in ('manual QA', ):
                continue
            if partype not in commands[sbname]:
                commands[sbname][partype] = replacements
            else:
                for spwsel, tcpars in replacements.items():
                    if spwsel in commands[sbname][partype]:
                        # we're replacing/overriding arguments here
                        for key, val in tcpars.items():
                            orig = commands[sbname][partype][spwsel][key]
                            if verbose:
                                print(f"{sbname} {partype} {spwsel}: {key}: {orig} -> {val}")
                            commands[sbname][partype][spwsel][key] = val
                    else:
                        # but if a spw was totally skipped, we replace it
                        # with the override version
                        if verbose:
                            print(f"SPW {spwsel} was skipped in {sbname} and is being replaced")
                        commands[sbname][partype][spwsel] = tcpars
    return commands


def main():

    with open(f"{pipedir}/default_tclean_commands.json", "r") as fh:
        default_commands = json.load(fh)

    commands = merge_aggregate(default_commands)
    commands = merge_override(commands)

    return commands


def get_commands():
    return main()


def make_table():
    clean_keys = list(set(itertools.chain.from_iterable(
        [list(commands[key]['tclean_cube_pars'][spw].keys())
         for key in commands
         for spw in commands[key]['tclean_cube_pars']
        ])))
    tabledata = [[key, spw,] + [(commands[key]['tclean_cube_pars'][spw][tkey]
                            if tkey in commands[key]['tclean_cube_pars'][spw]
                            else None)
                           for tkey in clean_keys]
                   for key in commands
                   for spw in commands[key]['tclean_cube_pars']
                  ]
    tabledata = list(map(list, zip(*tabledata)))
    table = Table(tabledata,
                  names=['fieldname', 'spwname',] + clean_keys )
    return table


if __name__ == "__main__":
    commands = main()
