"""
Environmental variables:

    ACES_ROOTDIR: the reduction_ACES directory
    ACES_DATADIR: the root data directory (e.g. /orange/adamginsburg/ACES/rawdata/)
    RUNONCE: boolean, if set, will only do one imaging run before quitting

Optional:
    PROJCODE:
    SOUS
    GOUS
"""
import os, sys, glob, json
# from casatasks import tclean
def tclean(**kwargs):
    """fake"""
    print(kwargs['imagename'], kwargs['parallel'])

if os.getenv('ACES_ROOTDIR') is None:
    raise ValueError("Specify ACES_ROOTDIR environment variable ")
else:
    rootdir = os.environ['ACES_ROOTDIR']
    sys.path.append(rootdir)
    sys.path.append(f'{rootdir}/pipeline_scripts')

# if this isn't in the env pars, we get an intentional crash:
# you have to specify that.
datadir = os.environ['ACES_DATADIR']

projcode = os.getenv('PROJCODE') or '2021.1.00172.L'
sous = os.getenv('SOUS') or 'A001_X1590_X30a8'
gous = os.getenv('GOUS') or 'A001_X1590_X30a9'
runonce = bool(os.getenv('RUNONCE'))

from merge_tclean_commands import commands

suffixes = {"tclean_cont_pars": ("image.tt0", "residual.tt0", "model.tt0", "psf.tt0"),
            "tclean_cube_pars": ("image", "residual", "model", "psf"),
            }

for sbname,allpars in commands.items():
    mous_ = allpars['mous']
    mous = mous_[6:].replace("/", "_")
    assert len(mous) in (14, 15)
    workingpath = f'{datadir}/{projcode}/science_goal.uid___{sous}/group.uid___{gous}/member.uid___{mous}/calibrated/working'
    if os.path.exists(workingpath):
        for partype in suffixes.keys():
            parsset = allpars[partype]
            for spwsel, tcpars in parsset.items():
                baseimname = tcpars['imagename']
                stage_wildcard_name = ".".join(baseimname.split(".")[:1] + ["*"] + baseimname.split(".")[2:])
                imname = f"{workingpath}/{baseimname}"
                exists = {suffix: os.path.exists(f"{imname}.{suffix}") for suffix in suffixes[partype]}
                imname_wild = f"{workingpath}/{stage_wildcard_name}"
                exists_wild = {suffix: len(glob.glob(f"{imname_wild}.{suffix}")) > 0 for suffix in suffixes[partype]}
                if not any(exists.values()) and not any(exists_wild.values()):
                    os.chdir(workingpath)
                    print(f"{sbname} {partype} {spwsel}: ", end=" ")
                    tclean(**tcpars)
                    if runonce:
                        sys.exit(0)
                elif all(exists_wild.values()):
                    #print(f"Found all files exist for {mous} but from a different stage")
                    pass # this is the "OK" state
                elif all(exists.values()):
                    #print(f"Found all files exist for {mous}")
                    pass # this is the "OK" state
                else:
                    print(f"Found partially-completed run for {sbname} {partype}: {exists_wild}")
    else:
        print(f"Did not find mous {mous}")

