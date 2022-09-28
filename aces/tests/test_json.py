import json
from aces import conf
rootdir = conf.basepath

with open(f"{rootdir}/reduction_ACES/aces/pipeline_scripts/default_tclean_commands.json", "r") as fh:
    default_commands = json.load(fh)

with open(f"{rootdir}/reduction_ACES/aces/pipeline_scripts/aggregate_high_tclean_commands.json", "r") as fh:
    aggregate_high_commands = json.load(fh)

with open(f"{rootdir}/reduction_ACES/aces/pipeline_scripts/override_tclean_commands.json", "r") as fh:
    override_commands = json.load(fh)
