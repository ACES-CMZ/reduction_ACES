import json
import pytest
from aces import conf
import os

# the tests have to run on tox
# rootdir = os.getcwd()
#rootdir = conf.basepath

rootdir = os.path.join(os.path.dirname(__file__), '..', 'pipeline_scripts')

with open(f"{rootdir}/default_tclean_commands.json", "r") as fh:
    default_commands = json.load(fh)

with open(f"{rootdir}/aggregate_high_tclean_commands.json", "r") as fh:
    aggregate_high_commands = json.load(fh)

with open(f"{rootdir}/override_tclean_commands.json", "r") as fh:
    override_commands = json.load(fh)


@pytest.mark.parametrize('jsonfile', ('default_tclean_commands.json',
                                      'aggregate_high_tclean_commands.json',
                                      'override_tclean_commands.json'))
def test_json(jsonfile):
    with open(f"{rootdir}/{jsonfile}", "r") as fh:
        override_commands = json.load(fh)
    assert len(override_commands) > 0
