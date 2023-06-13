import os
import textwrap
import datetime
from astropy import units as u
import subprocess

def parallel_clean_slurm(nchan, imagename, spw, start=0, width=1, nchan_per=128,
                         ntasks=4, mem_per_cpu='4gb', jobname='array_clean',
                         account='astronomy-dept', qos='astronomy-dept-b',
                         jobtime='96:00:00',
                         CASAVERSION = 'casa-6.4.3-2-pipeline-2021.3.0.17', #'casa-6.5.5-21-py3.8',
                         field='Sgr_A_star',
                         workdir='/blue/adamginsburg/adamginsburg/ACES/workdir',
                         dry=False,
                         savedir=None,
                         **kwargs):

    try:
        hasunit = False
        start = int(start)
        width = int(width)
    except ValueError:
        hasunit = True
        start = u.Quantity(start).to(u.GHz).value
        width = u.Quantity(width).to(u.GHz).value

    NARRAY = nchan // nchan_per
    assert NARRAY > 1

    tclean_kwargs = {'nchan': nchan_per,}
    tclean_kwargs.update(**kwargs)
    # can't be zero, even though docs say it can be
    del tclean_kwargs['interactive']
    del tclean_kwargs['parallel']
    assert 'interactive' not in tclean_kwargs
    tclean_kwargs['calcres'] = True
    tclean_kwargs['calcpsf'] = True


    # TODO:
    # since we're no longer splitting out subsections of the MS, the split should only be done once
    # and the split job should be a dependency of the tclean array job.
    # In the current version, there is probably a race condition.

    rename_vis = textwrap.dedent(f"""

        def rename_vis(x):
            return os.path.join("{workdir}",
                                os.path.basename(x).replace('.ms',
                                                            f'_spw{spw}.ms'))

                                 """)
    logprint = textwrap.dedent("""

        def logprint(x):
            print(x, flush=True)
            casalog.post(str(x), origin='parallel_tclean')

                               """)

    splitcmd = textwrap.dedent(
        f"""

        # to be called AFTER tclean_kwargs are set

        import sys

        def test_valid(vis):
            try:
                msmd.open(vis)
                nchan_max = msmd.nchan(0)
                msmd.close()
                return True
            except RuntimeError as ex:
                logprint(f"vis {{vis}} was invalid.  Removing.")
                assert 'orange' not in vis
                shutil.rmtree(vis)
                return False


        for vis in {tclean_kwargs['vis']}:
            outputvis=f'{{rename_vis(vis)}}'
            if not os.path.exists(outputvis) or not test_valid(vis):
                try:
                    logprint(f"Splitting {{vis}} with defaults")
                    split(vis=vis,
                        outputvis=outputvis,
                        field='{field}',
                        spw=splitspw)
                    if not os.path.exists(outputvis):
                        raise ValueError("Did not split")
                    else:
                        logprint(f"Splitting {{vis}} with default (CORRECTED) was successful")
                except Exception as ex:
                    logprint(f"Failed first attempt with exception {{ex}}")
                    logprint(f"Splitting {{vis}} with datacolumn='data'")
                    split(vis=vis,
                          outputvis=outputvis,
                          field='{field}',
                          datacolumn='data',
                          spw=splitspw)

        for vis in {tclean_kwargs['vis']}:
            outputvis=f'{{rename_vis(vis)}}'
            if not os.path.exists(outputvis):
                # fail!
                sys.exit(1)

            # test that vis is valid
            msmd.open(vis)
            nchan_max = msmd.nchan(0)
            msmd.close()

        """)

    script = textwrap.dedent(
        f"""
        import os, shutil
        os.chdir('{workdir}')
        tclean_kwargs = {tclean_kwargs}
        width = {width}
        nchan_per = {nchan_per}
        splitspw = {spw}

        """)
    script += logprint
    script += "\nlogprint(f'Log file name is {os.getenv(\"LOGFILENAME\")}')\n"
    script += rename_vis
    splitscript = script + splitcmd

    if hasunit:
        script += textwrap.dedent(f"""
            start = int(os.getenv('SLURM_ARRAY_TASK_ID')) * {nchan_per} * {width} + {start}
            tclean_kwargs['start'] = f'{{start}}GHz'
            tclean_kwargs['width'] = f'{{width}}GHz'
            startchan = int(os.getenv('SLURM_ARRAY_TASK_ID')) * {nchan_per}
            tclean_kwargs['imagename'] = os.path.basename(f"{imagename}.{{startchan:04d}}.{{nchan_per:03d}}")
            #splitspw = f'{spw}:{{start-width}}GHz~{{start+width*(nchan_per+1)}}GHz'
            """)
    else:
        script += textwrap.dedent(f"""
            startchan = start = int(os.getenv('SLURM_ARRAY_TASK_ID')) * {nchan_per} * {width} + {start}
            tclean_kwargs['start'] = start
            tclean_kwargs['width'] = width
            tclean_kwargs['imagename'] = os.path.basename(f"{imagename}.{{start:04d}}.{{nchan_per:03d}}")
            """)


    script += textwrap.dedent(f"""

        # rename the vises whether or not we image so cleanup will work
        tclean_kwargs['vis'] = [rename_vis(vis) for vis in tclean_kwargs['vis']]

        for vis in tclean_kwargs['vis']:
            msmd.open(vis)
            # assume spw=0
            nchan_max = msmd.nchan(0)
            msmd.close()

        if start >= nchan_max:
            logprint(f"Maximum number of channels is {{nchan_max}} and start={{start}}.  Quitting.")
            sys.exit(0)

        if os.path.exists(tclean_kwargs['imagename'] + ".image"):
            logprint("Already done with startchan={{startchan}}")
        else:
            logprint(f'tclean_kwargs: {{tclean_kwargs}}')
            logprint(tclean_kwargs['vis'])
            logprint(f"Cleaning with startchan={{startchan}}")

            tclean(**tclean_kwargs)

            if not os.path.exists(tclean_kwargs['imagename'] + ".image"):
                raise ValueError(f"FAILURE: image {{tclean_kwargs['imagename']}}.image was not produced")

        """)


    splitscriptname = os.path.join(workdir, f"{imagename}_parallel_split_script.py")
    with open(splitscriptname, 'w') as fh:
        fh.write(splitscript)

    now = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")

    runsplitcmd = ("#!/bin/bash\n"
                   f'LOGFILENAME="casa_log_split_{jobname}_${{SLURM_JOBID}}_${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}_{now}.log"\n'
                   'echo "Log file is $LOGFILENAME"\n'
                   f'xvfb-run -d /orange/adamginsburg/casa/{CASAVERSION}/bin/casa'
                   ' --nologger --nogui '
                   ' --logfile=${LOGFILENAME} '
                   f' -c "execfile(\'{splitscriptname}\')"')

    slurmsplitcmdsh = imagename+"_split_slurm_cmd.sh"
    with open(slurmsplitcmdsh, 'w') as fh:
        fh.write(runsplitcmd)

    slurmsplitcmd = (f'/opt/slurm/bin/sbatch --ntasks={ntasks} '
                     f'--mem-per-cpu={mem_per_cpu} --output={jobname}_%j_%A_%a.log '
                     f'--job-name={jobname}_split --account={account} '
                     f'--qos={qos} --export=ALL --time={jobtime} {slurmsplitcmdsh}\n')


    if dry:
        print(slurmsplitcmd.split())
        scriptjobid = 'PLACEHOLDER'
    else:
        print(slurmsplitcmd.split())
        sbatch = subprocess.check_output(slurmsplitcmd.split())
        scriptjobid = sbatch.decode().split()[-1]
        print(f'Split: {sbatch.decode()} with jobname={jobname}')


    scriptname = os.path.join(workdir, f"{imagename}_parallel_script.py")
    with open(scriptname, 'w') as fh:
        fh.write(script)



    runcmd = ("#!/bin/bash\n"
              f'LOGFILENAME="casa_log_line_{jobname}_${{SLURM_JOBID}}_${{SLURM_ARRAY_TASK_ID}}_{now}.log"\n'
              'echo "Log file is $LOGFILENAME"\n'
              f'xvfb-run -d /orange/adamginsburg/casa/{CASAVERSION}/bin/casa'
              ' --nologger --nogui '
              ' --logfile=${LOGFILENAME} '
              f' -c "execfile(\'{scriptname}\')"')

    slurmcmd = imagename+"_slurm_cmd.sh"
    with open(slurmcmd, 'w') as fh:
        fh.write(runcmd)

    cmd = (f'/opt/slurm/bin/sbatch --ntasks={ntasks} '
           f'--mem-per-cpu={mem_per_cpu} --output={jobname}_%j_%A_%a.log '
           f'--job-name={jobname}_arr --account={account} '
           f'--array=0-{NARRAY} '
           f'--dependency=afterok:{scriptjobid} '
           f'--qos={qos} --export=ALL --time={jobtime} {slurmcmd}\n')

    if dry:
        print(cmd.split())
        jobid = 'PLACEHOLDER'
    else:
        print(cmd.split())
        sbatch = subprocess.check_output(cmd.split())
        jobid = sbatch.decode().split()[-1]
        print(f'{sbatch.decode()} with jobname={jobname}')


    mergescriptname = os.path.join(workdir, imagename+"_merge_script.py")

    # note the forced dedenting here because of {logprint} and {rename_vis}
    mergescript = textwrap.dedent(
f"""
import glob, os, shutil, datetime
savedir = '{savedir}'
os.chdir('{workdir}')
for suffix in ("image", "pb", "psf", "model", "residual", "weight", "mask", "image.pbcor", "sumwt"):
    outfile = os.path.basename(f'{imagename}.{{suffix}}')
    infiles = sorted(glob.glob(os.path.basename(f'{imagename}.[0-9]*.{{suffix}}')))
    print(outfile, infiles)
    if suffix in ("image.pbcor", "sumwt"):
        # these may not always exist
        print(f"Found only {{len(infiles)}} files: {{infiles}}.  suffix={{suffix}}")
    else:
        assert len(infiles) > 0, f"Found only {{len(infiles)}} files: {{infiles}}.  suffix={{suffix}}"

for suffix in ("image", "pb", "psf", "model", "residual", "weight", "mask", "image.pbcor", "sumwt"):
    outfile = os.path.basename(f'{imagename}.{{suffix}}')
    infiles = sorted(glob.glob(os.path.basename(f'{imagename}.[0-9]*.{{suffix}}')))
    if len(infiles) == 0:
        print(f"Skipped suffix {{suffix}}")
        continue

    now = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
    if os.path.exists(outfile):
        backup_outfile = outfile+f".backup_{{now}}"
        print(f"Found existing file {{outfile}}.  Moving to {{new_outfile}}")
        shutil.move(outfile, backup_outfile)

    ia.imageconcat(outfile=outfile,
                   infiles=infiles,
                   mode='m')
    if savedir and os.path.exists(savedir):
        print(f"Moving {{outfile}} to {{savedir}}")
        full_outfile = os.path.join(savedir, outfile)
        if os.path.exists(full_outfile):
            print("Outfile {{full_outfile}} already exists.  Check what's up.")
        else:
            shutil.move(outfile, savedir)
    else:
        print(f"Savedir {{savedir}} does not exist")


# Cleanup stage

import os, shutil
os.chdir('{workdir}')
tclean_kwargs = {tclean_kwargs}

{logprint}
{rename_vis}

for vis in tclean_kwargs['vis']:
    vis = rename_vis(vis)
    logprint(f"Removing visibility {{vis}}")
    assert 'orange' not in vis
    shutil.rmtree(vis)
""")

    with open(mergescriptname, 'w') as fh:
        fh.write(mergescript)


    #now = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
    #LOGFILENAME = f"casa_log_line_{jobname}_merge_{now}.log"
    runcmd_merge = ("#!/bin/bash\n"
                    f'LOGFILENAME="casa_log_merge_{jobname}_${{SLURM_JOBID}}_{now}.log"\n'
                    'echo "Log file is $LOGFILENAME"\n'
                    f'xvfb-run -d /orange/adamginsburg/casa/{CASAVERSION}/bin/casa'
                    f' --nologger --nogui '
                    f' --logfile=$LOGFILENAME '
                    f' -c "execfile(\'{mergescriptname}\')"')

    slurmcmd_merge = imagename+"_slurm_cmd_merge.sh"
    with open(slurmcmd_merge, 'w') as fh:
        fh.write(runcmd_merge)

    cmd = (f'/opt/slurm/bin/sbatch --ntasks={ntasks*4} '
           f'--mem-per-cpu={mem_per_cpu} --output={jobname}_merge_%j_%A_%a.log --job-name={jobname}_merge --account={account} '
           f'--dependency=afterok:{jobid} '
           f'--qos={qos} --export=ALL --time={jobtime} {slurmcmd_merge}')

    if dry:
        print(cmd.split())
    else:
        print(cmd.split())
        sbatch = subprocess.check_output(cmd.split())
        print(f'{sbatch.decode()} with jobname={jobname}_merge')
