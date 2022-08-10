# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
# ----------------------------------------------------------------------------

__all__ = []


from astropy import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `aces`.
    """

    basepath = _config.ConfigItem(
        '/orange/adamginsburg/ACES/'
        'Base path in which data/ exists')


conf = Conf()

import retrieval_scripts

__all__ = ['retrieval_scripts',
           'Conf', 'conf',
           ]
