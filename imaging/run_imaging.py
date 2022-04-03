"""
Environmental variables:

    ACES_ROOTDIR: the reduction_ACES directory
    ACES_DATADIR: the root data directory (e.g. /orange/adamginsburg/ACES/rawdata/)
    RUNONCE: boolean, if set, will only do one imaging run before quitting
    DUMMYRUN: boolean, if set, will only print what it's going to do
    TRYDROPTARGET: boolean.  If set, will try MSes w/o "_target" in them if it can't find them with

Optional:
    PROJCODE:
    SOUS
    GOUS
"""
import os, sys, glob, json, shutil

if os.getenv('DUMMYRUN'):
    def tclean(**kwargs):
        """fake"""
        print(kwargs['imagename'], kwargs['parallel'])
else:
    from casatasks import tclean

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
cleanup = bool(os.getenv('CLEANUP'))
scriptlist = os.path.realpath(os.getenv('SCRIPTLIST') or './scriptlist.txt')
print(f"RUNONCE={runonce} CLEANUP={cleanup} DUMMYRUN={bool(os.getenv('DUMMYRUN'))} SCRIPTLIST={scriptlist}")

# touch scriptlist
with open(scriptlist, 'w') as fh:
    fh.write("")


from merge_tclean_commands import commands

suffixes = {"tclean_cont_pars": ("image.tt0", "residual.tt0", "model.tt0", "psf.tt0"),
            "tclean_cube_pars": ("image", "residual", "model", "psf"),
            }

# these aren't really user-configurable
tcpars_override = {'calcpsf': True, 'calcres': True,}

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
                    if not all(os.path.exists(x) for x in tcpars['vis']) and os.getenv('TRYDROPTARGET'):
                        tcpars['vis'] = [x.replace("_target", "") for x in tcpars["vis"]]
                    if not all(os.path.exists(x) for x in tcpars['vis']):
                        print(f"ERROR: Files not found: {tcpars['vis']}")
                        continue
                    tcpars.update(tcpars_override)

                    # try using full path?
                    tcpars['vis'] = [os.path.realpath(x) for x in tcpars["vis"]]
                    tcpars['imagename'] = os.path.realpath(tcpars['imagename'])

                    print(f"Creating script for {partype} tclean in {workingpath} for sb {sbname} with kwargs: \n{tcpars}")

                    with open(f"tclean_{partype}_{sbname}_{spwsel}.py", "w") as fh:
                        fh.write(f"tclean(\n")
                        for key, val in tcpars.items():
                            fh.write(f"       {key}={repr(val)},\n")
                        fh.write(")")

                    with open(scriptlist, 'a') as fh:
                        fh.write(f'{workingpath}/tclean_{partype}_{sbname}_{spwsel}.py\n')

                    #tclean(**tcpars)

                    #if runonce:
                    #    sys.exit(0)
                elif all(exists.values()):
                    #print(f"Found all files exist for {mous}")
                    pass # this is the "OK" state
                elif all(exists_wild.values()):
                    #print(f"Found all files exist for {mous} but from a different stage")
                    pass # this is the "OK" state
                else:
                    if any(exists.values()):
                        print(f"Found partially-completed run for {sbname} {partype} {spwsel}: {exists}")
                    elif any(exists_wild.values()):
                        print(f"Found partially-completed run for {sbname} {partype} {spwsel}: glob-based {exists_wild}")
                    else:
                        raise ValueError("Huh?  No.")
                    if cleanup:
                        for suffix in exists:
                            if exists[suffix]:
                                print(f"Removing existing file with suffix {suffix}: {imname}.{suffix}")
                                shutil.rmtree(f'{imname}.{suffix}')
                            elif exists_wild[suffix]:
                                print(f"DRY: Removing existing files with suffix {suffix}: ", glob.glob(f"{imname_wild}.{suffix}"))
    else:
        print(f"Did not find mous {mous}")

if runonce:
    print("Completed re-imaging run with RUNONCE enabled, but didn't run at all.")
