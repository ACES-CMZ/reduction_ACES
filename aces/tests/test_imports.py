# flake8: noqa
import aces
import aces.imaging

# from aces.joint_deconvolution import split_spws_7m
# from aces.joint_deconvolution import split_spws_12m
# from aces.joint_deconvolution import joint_deconvolution_and_feather_auto
# from aces.joint_deconvolution import scripts import parse_contdotdat
# from aces.joint_deconvolution import scripts import __init__
# from aces.joint_deconvolution import __init__

from aces.pipeline_scripts import generate_spw33_commands
from aces.pipeline_scripts import merge_tclean_commands
from aces.pipeline_scripts import recover_tclean_commands

# CASA, not valid python
#from aces.pipeline_scripts import imaging_pipeline_rerun

from aces.retrieval_scripts import parse_weblog
from aces.retrieval_scripts import make_humanreadable_links
from aces.retrieval_scripts import retrieve_weblogs
# CASA from aces.retrieval_scripts import run_pipeline
from aces.retrieval_scripts import retrieve_data
from aces.retrieval_scripts import mous_map

# CASA from aces.hipergator_scripts import hack_plotms
from aces.hipergator_scripts import link_repipeline_weblogs
from aces.hipergator_scripts import job_runner
from aces.hipergator_scripts import ghapi_update
from aces.hipergator_scripts import delivery_status

from aces.analysis import imstats
from aces.analysis import spectral_extraction_Feb2022
from aces.analysis import cube_stats_grid

from aces.observation_planning import spectral_shift_planning

from aces.imaging import write_tclean_scripts
from aces.imaging import make_mosaic
from aces.imaging import mosaic_7m
from aces.imaging import mosaic_12m
from aces.imaging import mosaic_TP


def test_imports():
    from aces import conf
    import os

    # the tests have to run on tox
    conf.basepath = os.getcwd()
    # set the directory to the parent of the parent of this one (so, should be the package-containing directory, under which aces/ will exist)
    conf.basepath = os.path.realpath(os.path.join(os.path.dirname(__file__), '../..'))
    print(conf.basepath)

    from aces.pipeline_scripts import generate_spw33_commands
    from aces.pipeline_scripts import merge_tclean_commands
    from aces.pipeline_scripts import recover_tclean_commands
    # CASA, not valid python
    #from aces.pipeline_scripts import imaging_pipeline_rerun

    from aces.retrieval_scripts import parse_weblog
    from aces.retrieval_scripts import make_humanreadable_links
    from aces.retrieval_scripts import retrieve_weblogs
    # CASA from aces.retrieval_scripts import run_pipeline
    from aces.retrieval_scripts import retrieve_data
    from aces.retrieval_scripts import mous_map

    # CASA from aces.hipergator_scripts import hack_plotms
    # don't test hipergator scripts unless
    import socket
    if 'ufhpc' in socket.gethostname():
        from aces.hipergator_scripts import link_repipeline_weblogs
        from aces.hipergator_scripts import job_runner
        from aces.hipergator_scripts import ghapi_update
        from aces.hipergator_scripts import delivery_status

    from aces.analysis import imstats
    from aces.analysis import spectral_extraction_Feb2022
    from aces.analysis import cube_stats_grid

    from aces.observation_planning import spectral_shift_planning

    from aces.imaging import write_tclean_scripts
    from aces.imaging import make_mosaic
    from aces.imaging import mosaic_7m
    from aces.imaging import mosaic_12m
    from aces.imaging import mosaic_TP