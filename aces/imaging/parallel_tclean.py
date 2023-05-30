import os
import textwrap
import datetime
from astropy import units as u
import subprocess

def parallel_clean_slurm(nchan, imagename, spw, start=0, width=1, nchan_per=128,
                         ntasks=4, mem_per_cpu='4gb', jobname='array_clean',
                         account='astronomy-dept', qos='astronomy-dept-b',
                         jobtime='96:00:00',
                         CASAVERSION = 'casa-6.5.5-21-py3.8',
                         field='Sgr_A_star',
                         workdir='/blue/adamginsburg/adamginsburg/ACES/workdir',
                         dry=False,
                         #savedir='/blue/adamginsburg/adamginsburg/ACES/workdir',
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

    splitcmd = textwrap.dedent(
        f"""
            # to be called AFTER tclean_kwargs are set
            def rename_vis(x):
                return os.path.join("{workdir}",
                                    os.path.basename(x).replace('.ms',
                                                                f'_spw{spw}_ch{{startchan}}+{{nchan_per}}.ms'))

            for vis in {tclean_kwargs['vis']}:
                outputvis=f'{{rename_vis(vis)}}'
                if not os.path.exists(outputvis):
                    try:
                        logprint(f"Splitting {{vis}} with defaults")
                        split(vis=vis,
                            outputvis=outputvis,
                            field='{field}',
                            spw={{splitspw}})
                        if not os.path.exists(outputvis):
                            raise ValueError("Did not split")
                        else:
                            logprint(f"Splitting {{vis}} with default (CORRECTED) was successful")
                    except Exception as ex:
                        logprint(ex)
                        logprint(f"Splitting {{vis}} with datacolumn='data'")
                        split(vis=vis,
                              outputvis=outputvis,
                              field='{field}',
                              datacolumn='data',
                              spw={{splitspw}})\n
                              """)

    script = textwrap.dedent(
        f"""
        import os
        os.chdir('{workdir}')
        tclean_kwargs = {kwargs}
    """)

    if hasunit:
        script += textwrap.dedent(f"""
        start = int(os.getenv('SLURM_ARRAY_TASK_ID')) * {nchan_per} * {width} + {start}
        tclean_kwargs['start'] = f'{{start}}GHz'
        tclean_kwargs['width'] = f'{width}GHz'
        startchan = int(os.getenv('SLURM_ARRAY_TASK_ID')) * {nchan_per}
        tclean_kwargs['imagename'] = os.path.basename(f"{imagename}.{{startchan:04d}}.{{nchan_per:03d}}")
        splitspw = f'{spw}:{{start-width}}GHz~{{start+width*(nchan_per+1)}}GHz'
        """)
    else:
        script += textwrap.dedent(f"""
        startchan = start = int(os.getenv('SLURM_ARRAY_TASK_ID')) * {nchan_per} * {width} + {start}
        tclean_kwargs['start'] = start
        tclean_kwargs['width'] = width
        tclean_kwargs['imagename'] = os.path.basename(f"{imagename}.{{start:04d}}.{{nchan_per:03d}}")
        splitspw = spw
        """)

    script += splitcmd

    script += textwrap.dedent(f"""
    tclean_kwargs['vis'] = [rename_vis(vis) for vis in tclean_kwargs['vis']]
    print(tclean_kwargs['vis'])

    tclean(**tclean_kwargs)\n""")

    scriptname = os.path.join(workdir, imagename+"_parallel_script.py")
    with open(scriptname, 'w') as fh:
        fh.write(script)

    now = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
    LOGFILENAME = f"casa_log_line_{jobname}_{now}.log"

    runcmd = ("#!/bin/bash\n"
              f'xvfb-run -d /orange/adamginsburg/casa/{CASAVERSION}/bin/casa'
              f' --nologger --nogui '
              f' --logfile={LOGFILENAME} '
              f' -c "execfile(\'{scriptname}\')"')

    slurmcmd = imagename+"_slurm_cmd.sh"
    with open(slurmcmd, 'w') as fh:
        fh.write(runcmd)

    cmd = (f'/opt/slurm/bin/sbatch --ntasks={ntasks} '
           f'--mem-per-cpu={mem_per_cpu} --output={jobname}_%j_%A_%a.log --job-name={jobname} --account={account} '
           f'--array=0-{NARRAY} '
           f'--qos={qos} --export=ALL --time={jobtime} {slurmcmd}')

    if dry:
        print(cmd.split())
        jobid = 'PLACEHOLDER'
    else:
        print(cmd.split())
        sbatch = subprocess.check_output(cmd.split())
        jobid = sbatch.decode().split()[-1]
        print(f'{sbatch.decode()} with jobname={jobname}')


    mergescriptname = os.path.join(workdir, imagename+"_merge_script.py")
    mergescript = textwrap.dedent(
        f"""
        import glob, os
        os.chdir('{workdir}')
        ia.imageconcat(outfile=os.path.basename('{imagename}.image',)
                       infiles=sorted(glob.glob(os.path.basename('{imagename}.[0-9]*.image'))),
                       mode='c')
        """)

    with open(mergescriptname, 'w') as fh:
        fh.write(mergescript)


    now = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
    LOGFILENAME = f"casa_log_line_{jobname}_merge_{now}.log"
    runcmd_merge = ("#!/bin/bash\n"
                    f'xvfb-run -d /orange/adamginsburg/casa/{CASAVERSION}/bin/casa'
                    f' --nologger --nogui '
                    f' --logfile={LOGFILENAME} '
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
