import json
from aces.pipeline_scripts import generate_aggregate_high_commands
from aces import conf

generate_aggregate_high_commands.main()

# old version
# if os.getenv('ACES_ROOTDIR') is None:
#     raise ValueError("Specify ACES_ROOTDIR environment variable ")
# else:
#     rootdir = os.environ['ACES_ROOTDIR']

if 'verbose' not in locals():
    verbose = False

pipedir = os.path.dirname(__file__)

with open(f"{pipedir}/default_tclean_commands.json", "r") as fh:
    default_commands = json.load(fh)

with open(f"{pipedir}/aggregate_high_tclean_commands.json", "r") as fh:
    aggregate_high_commands = json.load(fh)

with open(f"{pipedir}/override_tclean_commands.json", "r") as fh:
    override_commands = json.load(fh)

commands = default_commands


for sbname, allpars in aggregate_high_commands.items():
    for partype, replacements in allpars.items():
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


for sbname, allpars in override_commands.items():
    for partype, replacements in allpars.items():
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
