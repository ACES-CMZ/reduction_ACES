"""
Environmental variables:

    ACES_ROOTDIR: the reduction_ACES directory
    ACES_DATADIR: the root data directory (e.g. /orange/adamginsburg/ACES/rawdata/)
    RUNONCE: boolean, if set, will only do one imaging run before quitting
     (deprecated - not used any more)
    DUMMYRUN: boolean, if set, will only print what it's going to do
    TRYDROPTARGET: boolean.  If set, will try MSes w/o "_target" in them if it can't find them with

Optional:
    PROJCODE:
    SOUS
    GOUS
"""
import os
import glob
import shutil
import textwrap
import string
from aces import conf
from aces.pipeline_scripts.merge_tclean_commands import commands

if os.getenv('DUMMYRUN'):
    def tclean(**kwargs):
        """fake"""
        print(kwargs['imagename'], kwargs['parallel'])
else:
    pass
    # from casatasks import tclean


def main():
    import time
    t0 = time.time()
    print("Starting aces_write_tclean_scripts")
    # if this isn't in the env pars, we get an intentional crash:
    # you have to specify that.
    datadir = f'{conf.basepath}/data/'  # os.environ['ACES_DATADIR']

    projcode = os.getenv('PROJCODE') or '2021.1.00172.L'
    sous = os.getenv('SOUS') or 'A001_X1590_X30a8'
    gous = os.getenv('GOUS') or 'A001_X1590_X30a9'
    runonce = bool(os.getenv('RUNONCE'))
    cleanup = bool(os.getenv('CLEANUP'))
    scriptlist = os.path.realpath(os.getenv('SCRIPTLIST') or './scriptlist.txt')

    if os.getenv('TEMPORARY_WORKING_DIRECTORY'):
        temp_workdir = os.getenv('TEMPORARY_WORKING_DIRECTORY')
    else:
        temp_workdir = "/blue/adamginsburg/adamginsburg/ACES/workdir"
        # temp_workdir = False

    print(f"RUNONCE={runonce} CLEANUP={cleanup} DUMMYRUN={bool(os.getenv('DUMMYRUN'))} SCRIPTLIST={scriptlist} TEMP_WORKDIR={temp_workdir}")

    # touch scriptlist
    with open(scriptlist, 'w') as fh:
        fh.write("")

    suffixes = {"tclean_cont_pars": ("image.tt0", "residual.tt0", "model.tt0", "psf.tt0"),
                "tclean_cube_pars": ("image", "residual", "model", "psf"),
                }

    # these aren't really user-configurable
    tcpars_override = {'calcpsf': True, 'interactive': 0 }

    for sbname, allpars in commands.items():
        mous_ = allpars['mous']
        mous = mous_[6:].replace("/", "_")
        assert len(mous) in (14, 15, 16)
        workingpath = f'{datadir}/{projcode}/science_goal.uid___{sous}/group.uid___{gous}/member.uid___{mous}/calibrated/working'
        if os.path.exists(workingpath):
            for partype in suffixes.keys():
                parsset = allpars[partype]
                for spwsel, tcpars in parsset.items():
                    if 'imagename' not in tcpars:
                        print(f"***** BROKEN  spw {spwsel} *****")
                        continue
                    baseimname = tcpars['imagename']
                    stage_wildcard_name = ".".join(baseimname.split(".")[:1] + ["*"] + baseimname.split(".")[2:])
                    imname = f"{workingpath}/{baseimname}"
                    exists = {suffix: os.path.exists(f"{imname}.{suffix}") for suffix in suffixes[partype]}
                    imname_wild = f"{workingpath}/{stage_wildcard_name}"
                    exists_wild = {suffix: len(glob.glob(f"{imname_wild}.{suffix}")) > 0 for suffix in suffixes[partype]}

                    # make the clean scripts
                    os.chdir(workingpath)

                    field = sbname.split("_")[3]
                    config = sbname.replace("_updated", "").split("_")[5]
                    print(f"{os.getcwd()} {sbname} {mous} {partype} {spwsel} {field} {config}: ", end=" ")
                    assert config in ('7M', 'TM1', 'TP')
                    if not all(os.path.exists(x) for x in tcpars['vis']) and os.getenv('TRYDROPTARGET'):
                        tcpars['vis'] = [x.replace("_targets", "") for x in tcpars["vis"]]
                        tcpars['vis'] = [x.replace("_target", "") for x in tcpars["vis"]]
                        tcpars['datacolumn'] = 'corrected'

                    # Nov 10, 2022: try removing "_lines" from files:
                    if all(os.path.exists(x.replace("_targets_line", "")) for x in tcpars['vis']):
                        tcpars['vis'] = [x.replace("_targets_line", "") for x in tcpars["vis"]]
                    if all(os.path.exists(x.replace("_line", "")) for x in tcpars['vis']):
                        tcpars['vis'] = [x.replace("_line", "") for x in tcpars["vis"]]

                    if not all(os.path.exists(x) for x in tcpars['vis']):
                        print(f"ERROR: Files not found: {tcpars['vis']}")
                        continue
                    tcpars.update(tcpars_override)

                    # use full path for vis
                    tcpars['vis'] = [os.path.realpath(x) for x in tcpars["vis"]]

                    # but use *relative* path for imagename so we can run this in temp directories
                    orig_dirname = os.path.dirname(tcpars['imagename'])
                    tcpars['imagename'] = os.path.basename(tcpars['imagename'])
                    print(f"tcpars['imagename'] = {tcpars['imagename']}")

                    # dirname is a bit redundant now
                    dirname = '.'

                    imtype = tcpars['specmode']

                    tempdir_name = f'{temp_workdir}/{field}_{spwsel}_{imtype}_{config}_{mous}'

                    if not os.path.exists(tempdir_name) or not os.path.isdir(tempdir_name):
                        os.mkdir(tempdir_name)

                    # save directory is same as tempdir unless SLURM_TMPDIR is used
                    savedir_name = tempdir_name

                    # check for PSF
                    expected_psfname = os.path.join(tempdir_name if temp_workdir else dirname,
                                                    os.path.basename(tcpars['imagename']) +
                                                    ('.psf.tt0' if tcpars['specmode'] == 'mfs' else '.psf')
                                                    )
                    check_psf_exists = textwrap.dedent(f"""
                                            # if the PSF exists, don't re-calculate it
                                            calcpsf = not os.path.exists('{expected_psfname}')

                                            if not os.path.exists('{os.path.basename(expected_psfname)}') and not calcpsf:
                                                shutil.copytree('{expected_psfname}', '{os.path.basename(expected_psfname)}')
                                            \n\n""")
                    # this is overridden below tcpars['calcpsf'] = 'calcpsf'

                    if temp_workdir:

                        if 'aggregate' in spwsel:
                            def rename(x):
                                return os.path.join(tempdir_name, os.path.basename(x))

                            splitcmd = copycmds = "\n".join(
                                ["import shutil"] +
                                [textwrap.dedent(
                                    f"""
                                    try:
                                        shutil.copytree('{x}', '{rename(x)}')
                                    except FileExistsError as ex:
                                        print(f'MS file already copied: {{ex}}.  Proceeding.')
                                    """)
                                 for x in tcpars['vis']
                                 ])

                            # hard code that parallel = False for non-MPI runs
                            tcpars['parallel'] = False
                        else:
                            spw = int(spwsel.lstrip('spw'))

                            def rename(x):
                                return os.path.join(tempdir_name,
                                                    os.path.basename(x).replace('.ms', f'_spw{spw}.ms'))
                            splitcmd = [textwrap.dedent(
                                f"""
                                    outputvis='{rename(vis)}'
                                    if not os.path.exists(outputvis):
                                        try:
                                            logprint("Splitting {vis} with defaults")
                                            split(vis='{vis}',
                                                outputvis=outputvis,
                                                field='Sgr_A_star',
                                                spw={spw})
                                            if not os.path.exists(outputvis):
                                                raise ValueError("Did not split")
                                            else:
                                                logprint("Splitting {vis} with default (CORRECTED) was successful")
                                        except Exception as ex:
                                            logprint(ex)
                                            logprint("Splitting {vis} with datacolumn='data'")
                                            split(vis='{vis}',
                                                outputvis=outputvis,
                                                field='Sgr_A_star',
                                                datacolumn='data',
                                                spw={spw})
                                                """) for vis in tcpars['vis']]

                            # ONLY for line cubes, which are individually split out:
                            # the spw selection should now be 'everything in the MS'
                            tcpars['spw'] = ''

                        # all 'vis' must be renamed because of their new locations
                        tcpars['vis'] = [rename(x) for x in tcpars["vis"]]

                        savecmds = textwrap.dedent(
                             f"""
                             import numpy as np
                             import glob
                             def savedata():
                                 flist = glob.glob('{os.path.basename(tcpars['imagename'])}.*')
                                 for fn in flist:
                                     logprint(f'Copying {{fn}} to {savedir_name}')
                                     target = f'{savedir_name}/{{os.path.basename(fn)}}'
                                     if os.path.realpath(fn) == os.path.realpath(target):
                                        print("Skipping copy - source = destination")
                                     else:
                                         if fn.endswith('.fits'):
                                             shutil.copy(fn, target)
                                         else:
                                             if os.path.exists(target):
                                                 logprint(f'Removing {{target}} because it exists')
                                                 assert 'iter1' in f'{{os.path.basename(fn)}}'  # sanity check - don't remove important directories!
                                                 shutil.rmtree(target)
                                             shutil.copytree(fn, target)\n\n
                             """)

                        cleanupcmds = (textwrap.dedent(
                             f"""
                             import glob
                             flist = glob.glob('{workingpath}/{os.path.basename(tcpars['imagename'])}.*')
                             for fn in flist:
                                 logprint(f'Moving {{fn}} to {workingpath}')
                                 target = f'{workingpath}/{{os.path.basename(fn)}}'
                                 if os.path.exists(target):
                                     logprint(f'Removing {{target}} because it exists')
                                     assert 'iter1' in target  # sanity check - don't remove important directories!
                                     if fn.endswith('.fits'):
                                         os.remove(target)
                                     else:
                                         shutil.rmtree(target)
                                 shutil.move(fn, '{workingpath}/')\n\n""") +
                             "\n".join([f"shutil.rmtree('{tempdir_name}/{os.path.basename(x)}')" for x in tcpars['vis']])
                        )
                        # use local name instead
                        #tcpars['imagename'] = os.path.join(tempdir_name, os.path.basename(tcpars['imagename']))

                    print(f"Creating script for {partype} tclean in {workingpath} for sb {sbname} ")
                    # with kwargs: \n{tcpars}")

                    with open(f"{partype}_{sbname}_{spwsel}.py", "w") as fh:
                        fh.write("import os, shutil, glob, string\n")
                        fh.write(textwrap.dedent(f"""
                                 try:
                                     from taskinit import casalog
                                 except ImportError:
                                     from casatasks import casalog

                                 def logprint(string):
                                     casalog.post(string, origin='tclean_script')
                                     print(string)

                                 mpi_ntasks = os.getenv('mpi_ntasks')
                                 if mpi_ntasks is not None:
                                     parallel = int(mpi_ntasks) > 1
                                 else:
                                     parallel = False

                                 tempdir_name = os.getenv("SLURM_TMPDIR")
                                 if tempdir_name is None:
                                     # this is grabbed from write_tclean_scripts
                                     tempdir_name = "{tempdir_name}"
                                 savedir_name = '{savedir_name}'

                                 if not os.path.exists(tempdir_name):
                                     os.mkdir(tempdir_name)
                                 os.chdir(tempdir_name)

                                 logprint(f"Temporary directory used is {{tempdir_name}}")

                                 """))
                        fh.write("logprint(f'Started CASA in {os.getcwd()}')\n")

                        # tclean
                        if temp_workdir:
                            fh.write("".join(splitcmd))
                            fh.write("\n\n")
                        fh.write(check_psf_exists)
                        print(f"tcpars['imagename'] = {tcpars['imagename']}")
                        fh.write(savecmds)
                        fh.write("tclean_pars = dict(\n")
                        for key, val in tcpars.items():
                            if key in ('parallel', ):
                                fh.write(f"       {key}={key},\n")
                            elif key not in ('calcpsf', 'calcres', 'nmajor', 'interactive'):
                                fh.write(f"       {key}={repr(val)},\n")
                            else:
                                if key in ('calcpsf', 'calcres', 'nmajor', 'interactive'):
                                    pass
                                else:
                                    raise ValueError(f"ERROR: encountered invalid / overridden tclean kwarg {key}:{val}")
                        fh.write(")\n\n\n")

                        threshold = float(tcpars['threshold'].strip(string.ascii_letters))

                        # first major cycle
                        fh.write(f"ret = tclean(nmajor=1, calcpsf=calcpsf, fullsummary=True, interactive=False, **tclean_pars)\n\n")
                        fh.write(textwrap.dedent("""
                                     peakres = 0
                                     for val1 in ret['summaryminor'].values():
                                         for val2 in val1.values():
                                             for val3 in val2.values():
                                                 peakres = max([peakres, np.max(val3['peakRes'])])
                                                 """))

                        # remaining major cycles
                        fh.write(textwrap.dedent(f"""
                                 nmajors = 1
                                 while peakres > {threshold}:
                                     peakres = 0
                                     for val1 in ret['summaryminor'].values():
                                         for val2 in val1.values():
                                             for val3 in val2.values():
                                                 peakres = max([peakres, val3['peakRes']])
                                     print(f"{{nmajors}}: Residual={{peakres}} > threshold {threshold}")
                                     nmajors += 1
                                     ret = tclean(nmajor=1,
                                                  calcpsf=False,
                                                  interactive=False,
                                                  fullsummary=True,
                                                  calcres=True, # sadly must always calcres, even when redundant
                                                  **tclean_pars)
                                     savedata()\n
                                 """)
                                )

                        expected_imname = (os.path.basename(tcpars['imagename']) +
                                           ('.image.tt0.pbcor'
                                            if tcpars['specmode'] == 'mfs'
                                            else '.image.pbcor')
                                          )

                        check_exists = textwrap.dedent(f"""
                                              if not os.path.exists('{expected_imname}'):
                                                  raise IOError('Expected output file {expected_imname} does not exist.')
                                                  sys.exit(1)
                                                  \n\n""")
                        fh.write(check_exists)

                        if tcpars['specmode'] == 'cube':
                            fh.write(f"exportfits('{tcpars['imagename']}.image.pbcor', '{tcpars['imagename']}.image.pbcor.fits', overwrite=True)\n\n\n")
                        elif tcpars['specmode'] == 'mfs':
                            fh.write(f"exportfits('{tcpars['imagename']}.image.tt0.pbcor', '{tcpars['imagename']}.image.tt0.pbcor.fits', overwrite=True)\n\n\n")
                        else:
                            raise ValueError(f"Specmode was neither cube nor mfs: specmode={tcpars['specmode']}")
                        if temp_workdir:
                            fh.write(cleanupcmds)

                    with open(scriptlist, 'a') as fh:
                        fh.write(f'{workingpath}/{partype}_{sbname}_{spwsel}.py\n')

                    if not all(exists.values()) and not all(exists_wild.values()):
                        pass
                        # tclean(**tcpars)

                        # if runonce:
                        #    sys.exit(0)
                    elif all(exists.values()):
                        # print(f"Found all files exist for {mous}")
                        pass  # this is the "OK" state
                    elif all(exists_wild.values()):
                        # print(f"Found all files exist for {mous} but from a different stage")
                        pass  # this is the "OK" state
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

                    # reset imagename back to original
                    tcpars['imagename'] = os.path.join(orig_dirname, os.path.basename(tcpars['imagename']))
        else:
            print(f"Did not find mous {mous}")

    if runonce:
        print("Completed re-imaging run with RUNONCE enabled, but didn't run at all.")

    print(f"Done with aces_write_tclean_scripts after t={time.time()-t0}")
    globals().update(locals())


if __name__ == "__main__":
    main()
