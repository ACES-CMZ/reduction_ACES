# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
from . import retrieval_scripts
from . import pipeline_scripts
from . import imaging
# ----------------------------------------------------------------------------

from astropy import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `aces`.
    """

    basepath = _config.ConfigItem(
        '/orange/adamginsburg/ACES/',
        'Base path in which data/ exists')
    logpath = _config.ConfigItem(
        '/blue/adamginsburg/adamginsburg/ACES/logs/',
        'High-performance system where logs can be written')
    workpath = _config.ConfigItem(
        '/blue/adamginsburg/adamginsburg/ACES/workdir',
        'High-performance system where intermediate data products can go')


conf = Conf()


__all__ = ['retrieval_scripts',
           'pipeline_scripts',
           'imaging',
           'Conf', 'conf',
           ]
