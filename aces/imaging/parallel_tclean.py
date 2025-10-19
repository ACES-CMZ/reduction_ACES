import os
import textwrap
import datetime
from astropy import units as u
import subprocess


def parallel_clean_slurm(nchan, imagename, spw, start=0, width=1, nchan_per=128,
                         ntasks=4, mem_per_cpu='4gb', jobname='array_clean',
                         account='astronomy-dept', qos='astronomy-dept-b',
                         jobtime='96:00:00',
                         CASAVERSION='casa-6.4.3-2-pipeline-2021.3.0.17',
                         field='Sgr_A_star',
                         workdir='/red/adamginsburg/ACES/workdir',
                         logdir='/red/adamginsburg/ACES/logs',
                         dry=False,
                         savedir=None,
                         remove_incomplete_psf=True,
                         remove_incomplete_weight=True,
                         suffixes_to_merge_and_export=None,
                         array_jobs=None,
                         merge_only=False,
                         **kwargs):
    """
    Parameters
    ----------
    savedir:
        Where to put the files in the end
    workdir:
        Where to store the intermediate files
    """

    print(f"Starting parallel clean in workdir={workdir} with casa={CASAVERSION}", flush=True)

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
    if array_jobs is None:
        array_jobs = f'0-{NARRAY}'

    tclean_kwargs = {'nchan': nchan_per, }
    tclean_kwargs.update(**kwargs)
    # can't be zero, even though docs say it can be
    if 'interactive' in tclean_kwargs:
        del tclean_kwargs['interactive']
    if 'parallel' in tclean_kwargs:
        del tclean_kwargs['parallel']
    assert 'interactive' not in tclean_kwargs
    tclean_kwargs['calcres'] = True
    tclean_kwargs['calcpsf'] = True

    if tclean_kwargs.get('gridder') == 'mosaic':
        # is mask necessary?  don't think so?
        necessary_suffixes = ("image", "pb", "psf", "model", "residual", "weight")
    else:
        necessary_suffixes = ("image", "pb", "psf", "model", "residual")

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
            if not os.path.exists(outputvis) or not test_valid(outputvis):
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
            listobs(vis=vis, listfile=f'{{vis}}.listobs', overwrite=True)
            listobs(vis=outputvis, listfile=f'{{outputvis}}.listobs', overwrite=True)

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
            logprint(f"Already done with startchan={{startchan}}")
            sys.exit(0)
        elif os.path.exists(tclean_kwargs['imagename'] + ".residual"):
            logprint(ValueError(f"{{tclean_kwargs['imagename']}}.residual exists.  Current state unclear."))
            #sys.exit(0)
            # assumption is: residual got made, but no major cycles completed
            logprint("Attempting to continue anyway.")
        elif os.path.exists(tclean_kwargs['imagename'] + ".model"):
            # assumption is: could still be running
            logprint(ValueError(f"{{tclean_kwargs['imagename']}}.model exists.  Current state unclear."))
            logprint("Quitting.")
            sys.exit(0)
        elif os.path.exists(tclean_kwargs['imagename'] + ".psf"):
            if {remove_incomplete_psf}:
                shutil.rmtree(tclean_kwargs['imagename'] + ".psf")
            else:
                raise ValueError(f"{{tclean_kwargs['imagename']}}.psf exists.  Remove it before continuing.")
        elif os.path.exists(tclean_kwargs['imagename'] + ".weight"):
            if {remove_incomplete_weight}:
                shutil.rmtree(tclean_kwargs['imagename'] + ".weight")
            else:
                raise ValueError(f"{{tclean_kwargs['imagename']}}.weight exists.  Remove it before continuing.")

        # if we're continuing from a partially-completed run
        # we still want calcres=True in case model components were made
        if os.path.exists(tclean_kwargs['imagename'] + ".psf"):
            tclean_kwargs['calcpsf'] = False

        logprint(f'tclean_kwargs: {{tclean_kwargs}}')
        logprint(tclean_kwargs['vis'])
        logprint(f"Cleaning with startchan={{startchan}}")

        tclean(**tclean_kwargs)

        if not os.path.exists(tclean_kwargs['imagename'] + ".image"):
            raise ValueError(f"FAILURE: image {{tclean_kwargs['imagename']}}.image was not produced")

        def check_file(fn):
            " ensure that nonzero values were written, otherwise the imaging task failed "
            from spectral_cube import SpectralCube
            import numpy as np
            cube = SpectralCube.read(fn, format='fits' if fn.endswith('.fits') else 'casa_image', use_dask=True)
            mx = cube.max()
            if mx == 0 or np.isnan(mx):
                raise ValueError(f"File {{fn}} has a maximum value of {{mx}}")

        check_file(tclean_kwargs['imagename'] + ".image")
        check_file(tclean_kwargs['imagename'] + ".residual")

        """)

    splitscriptname = os.path.join(workdir, f"{imagename}_parallel_split_script.py")
    with open(splitscriptname, 'w') as fh:
        fh.write(splitscript)
    print(f"Wrote splitscript {splitscriptname}")

    now = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")

    runsplitcmd = ("#!/bin/bash\n"
                   f'LOGFILENAME="{logdir}/casa_log_split_{jobname}_${{SLURM_JOBID}}_${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}_{now}.log"\n'
                   'echo "Log file is $LOGFILENAME"\n'
                   f'xvfb-run -d /orange/adamginsburg/casa/{CASAVERSION}/bin/casa'
                   ' --nologger --nogui '
                   ' --logfile=${LOGFILENAME} '
                   f' -c "execfile(\'{splitscriptname}\')"')

    slurmsplitcmdsh = imagename + "_split_slurm_cmd.sh"
    with open(slurmsplitcmdsh, 'w') as fh:
        fh.write(runsplitcmd)
    print(f"Wrote runsplit {slurmsplitcmdsh}")

    slurmsplitcmd = (f'/opt/slurm/bin/sbatch --ntasks={ntasks} '
                     f'--mem-per-cpu={mem_per_cpu} --output={logdir}/{jobname}_%j_%A_%a.log '
                     f'--job-name={jobname}_split --account={account} '
                     f'--qos={qos} --export=ALL --time={jobtime} {slurmsplitcmdsh}\n')


    if dry:
        print(slurmsplitcmd.split())
        scriptjobid = 'PLACEHOLDER'
    elif merge_only:
        print("Skipping split job")
        scriptjobid = 'PLACEHOLDER'
    else:
        print(slurmsplitcmd.split())
        sbatch = subprocess.check_output(slurmsplitcmd.split())
        scriptjobid = sbatch.decode().split()[-1]
        print(f'Split: {sbatch.decode()} with jobname={jobname}')


    scriptname = os.path.join(workdir, f"{imagename}_parallel_script.py")
    with open(scriptname, 'w') as fh:
        fh.write(script)
    print(f"Wrote script {scriptname}")



    runcmd = ("#!/bin/bash\n"
              f'LOGFILENAME="{logdir}/casa_log_line_{jobname}_${{SLURM_JOBID}}_${{SLURM_ARRAY_TASK_ID}}_{now}.log"\n'
              'echo "Log file is $LOGFILENAME"\n'
              f'xvfb-run -d /orange/adamginsburg/casa/{CASAVERSION}/bin/casa'
              ' --nologger --nogui '
              ' --logfile=${LOGFILENAME} '
              f' -c "execfile(\'{scriptname}\')"')

    slurmcmd = imagename + "_slurm_cmd.sh"
    with open(slurmcmd, 'w') as fh:
        fh.write(runcmd)
    print(f"Wrote command {slurmcmd}")

    cmd = (f'/opt/slurm/bin/sbatch --ntasks={ntasks} '
           f'--mem-per-cpu={mem_per_cpu} --output={logdir}/{jobname}_%j_%A_%a.log '
           f'--job-name={jobname}_arr --account={account} '
           f'--array={array_jobs} '
           f'--dependency=afterok:{scriptjobid} '
           f'--qos={qos} --export=ALL --time={jobtime} {slurmcmd}\n')

    if dry:
        print(cmd.split())
        jobid = 'PLACEHOLDER'
    elif merge_only:
        print("Skipping array jobs")
        jobid = 'PLACEHOLDER'
    else:
        print(cmd.split())
        sbatch = subprocess.check_output(cmd.split())
        jobid = sbatch.decode().split()[-1]
        print(f'{sbatch.decode()} with jobname={jobname}')


    mergescriptname = os.path.join(workdir, imagename + "_merge_script.py")

    # note the forced dedenting here because of {logprint} and {rename_vis}
    mergescript = textwrap.dedent(
f"""
import glob, os, shutil, datetime

{logprint}

suffixes_to_merge_and_export = {suffixes_to_merge_and_export}
savedir = '{savedir}'
os.chdir('{workdir}')

necessary_suffixes = {necessary_suffixes}
print(f"Necessary suffixes: {{necessary_suffixes}}", flush=True)

for suffix in necessary_suffixes:
    outfile = os.path.basename(f'{imagename}.{{suffix}}')
    infiles = sorted(glob.glob(os.path.basename(f'{imagename}.[0-9]*.{{suffix}}')))
    print(outfile, infiles)
    if suffix in ("image.pbcor", "sumwt"):
        # these may not always exist
        print(f"Found only {{len(infiles)}} files: {{infiles}}.  suffix={{suffix}}")
    else:
        assert len(infiles) > 0, f"Found only {{len(infiles)}} files: {{infiles}}.  suffix={{suffix}}"

if suffixes_to_merge_and_export is None and os.getenv('SUFFIXES_TO_MERGE_AND_EXPORT') is not None:
    suffixes_to_merge_and_export = os.getenv('SUFFIXES_TO_MERGE_AND_EXPORT').split(',')
elif suffixes_to_merge_and_export is None:
    suffixes_to_merge_and_export = {necessary_suffixes}

print(f"Suffixes to merge and export: {{suffixes_to_merge_and_export}}", flush=True)

for suffix in necessary_suffixes:
    outfile = os.path.basename(f'{imagename}.{{suffix}}')
    #infiles = sorted(glob.glob(os.path.basename(f'{imagename}.[0-9]*.{{suffix}}')))
    infiles = [f'{imagename}.{{start:04d}}.{nchan_per:03d}.{{suffix}}'
               for start in range(0, {nchan}, {nchan_per})]
    if len(infiles) == 0:
        print(f"Skipped suffix {{suffix}}")
        continue
    for fn in infiles:
        if not os.path.exists(fn):
            print(f"Failure: file {{fn}} did not exist")

    now = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
    if os.path.exists(outfile):
        backup_outfile = outfile+f".backup_{{now}}"
        print(f"Found existing file {{outfile}}.  Moving to {{backup_outfile}}")
        os.rename(outfile, backup_outfile)
    if os.path.exists(outfile+".move"):
        backup_outfile = outfile+".move"+f".backup_{{now}}"
        print(f"Found existing file {{outfile+'.move'}}.  Moving to {{backup_outfile}}")
        os.rename(outfile+".move", backup_outfile)

    # reads much more efficiently, is consistent with other cubes
    print(f"Concatenating length {{len(infiles)}} files {{infiles}} to {{outfile}}", flush=True)
    ia.imageconcat(outfile=outfile,
                   infiles=infiles,
                   mode='p')
    ia.close()


    # consolidates the remainder into a single virtual cube
    # (these files can be removed later)
    #ia.imageconcat(outfile=outfile+".move",
    #               infiles=infiles,
    #               mode='m')

    exportfits(imagename=outfile,
               fitsimage=outfile+".fits",
               overwrite=True # don't want to crash here, and don't expect FITS files to be hanging around...
               )

psffile = os.path.basename(f'{imagename}.psf')
ia.open(psffile)
commonbeam = ia.commonbeam()
ia.close()

imagefile = os.path.basename(f'{imagename}.image')
ia.open(imagefile)
rbeam = ia.restoringbeam()
ia.close()

# if any beam is not the common beam...
def check_manybeam(rbeam, commonbeam):
    if 'beams' in rbeam:
        for beam in rbeam['beams'].values():
            if beam['*0']['major'] != commonbeam['major'] or beam['*0']['minor'] != commonbeam['minor']:
                print(f"Manybeam: beam={{beam['*0']}}, commonbeam={{commonbeam}}")
                return True
    return False

manybeam = check_manybeam(rbeam, commonbeam)

if os.path.exists('{imagename}.image.multibeam') and os.path.exists('{imagename}.image'):
    print(f"Found {imagename}.image.multibeam and {imagename}.image.")
    if manybeam:
        print("However, the image has multiple beams, so it is not a commonbeam image.")
        print("Removing {imagename}.image.multibeam so that .image can be moved to .image.multibeam")
        shutil.rmtree('{imagename}.image.multibeam')

def imsmooth_spectral_cube(imagename, outfile, commonbeam):
    '''
    Use spectral-cube to imsmooth without loading the whole thing into memory
    '''

    from spectral_cube import SpectralCube
    import radio_beam
    from astropy import units as u

    # dask is used here; no choice
    cube = SpectralCube.read('{imagename}.model', format='casa_image')
    # set the beam size to the pixel scale
    try:
        pixscale = cube.wcs.proj_plane_pixel_scales()[0]
    except AttributeError:
        pixscale = np.abs(cube.wcs.wcs.cdelt[0]) * u.deg
    cube = cube.with_beam(radio_beam.Beam(major=pixscale))
    cbeam = radio_beam.Beam(major=commonbeam['major']['value']*u.Unit(commonbeam['major']['unit']),
                            minor=commonbeam['minor']['value']*u.Unit(commonbeam['minor']['unit']),
                            pa=commonbeam['pa']['value']*u.Unit(commonbeam['pa']['unit']))
    convcube = cube.convolve_to(cbeam)

    convcube.write(outfile + ".fits", format='fits', overwrite=True)
    importfits(fitsimage=outfile + ".fits", imagename=outfile)


if manybeam:
    print("Beginning manybeam", flush=True)
    try:
        os.rename('{imagename}.image', '{imagename}.image.multibeam')
        print(f"Successfully moved {imagename}.image -> {imagename}.image.multibeam")
    except Exception as ex:
        print("Failed to move {imagename}.image -> {imagename}.image.multibeam, probably because the latter exists")
        print(ex)
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
        os.rename('{imagename}.image', f'{imagename}.image.{{timestamp}}')
        print(f"INSTEAD moved {imagename}.image -> {imagename}.image.{{timestamp}}")
        # the above is needed because we do _not_ want to declare victory when there are multibeam images

    assert not os.path.exists('{imagename}.image'), "FAILURE: {imagename}.image exists after moving it to {imagename}.image.multibeam"

    # it's possible for pbcor to not exist if we're not doing pbcor (to save space)
    if os.path.exists('{imagename}.image.pbcor'):
        try:
            os.rename('{imagename}.image.pbcor', '{imagename}.image.pbcor.multibeam')
            print(f"Successfully moved {imagename}.image.pbcor -> {imagename}.image.pbcor.multibeam")
        except Exception as ex:
            print("Failed to move {imagename}.image.pbcor -> {imagename}.image.pbcor.multibeam, probably because the latter exists")
            print(ex)
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
            os.rename('{imagename}.image.pbcor', f'{imagename}.image.pbcor.{{timestamp}}')
            print(f"INSTEAD moved {imagename}.image.pbcor -> {imagename}.image.pbcor.{{timestamp}}")

    if not os.path.exists('{imagename}.convmodel'):
        # Convmodel frequenty runs out of memory, so we use spectral-cube instead
        logprint('Creating convmodel {os.path.basename(imagename)}.convmodel')
        imsmooth_spectral_cube(imagename='{imagename}.model',
                               outfile='{imagename}.convmodel',
                               commonbeam=commonbeam)
        print(f"Successfully created convolved model {os.path.basename(imagename)}.convmodel: {{os.path.exists('{imagename}.convmodel')}}")
    if not os.path.exists('{os.path.basename(imagename)}.image'):
        logprint('Creating image {os.path.basename(imagename)}.image (it did not exist)')
        try:
            ia.imagecalc(outfile='{os.path.basename(imagename)}.image',
                        pixels='"{os.path.basename(imagename)}.convmodel" + "{os.path.basename(imagename)}.residual"',
                        imagemd='{os.path.basename(imagename)}.convmodel',
                        overwrite=False)
            ia.close()
        except RuntimeError as ex:
            ia.close()
            print(f"RuntimeError: {{ex}} - trying again with overwrite=False and to .image-temp instead")
            ia.imagecalc(outfile='{os.path.basename(imagename)}.image-temp',
                        pixels='"{os.path.basename(imagename)}.convmodel" + "{os.path.basename(imagename)}.residual"',
                        imagemd='{os.path.basename(imagename)}.convmodel',
                        overwrite=False)
            ia.close()
            if os.path.exists('{imagename}.image'):
                shutil.rmtree('{imagename}.image')
            shutil.move('{imagename}.image-temp', '{imagename}.image')
        print(f"Successfully created {os.path.basename(imagename)}.image: {{os.path.exists('{imagename}.image')}}")
    if not os.path.exists('{os.path.basename(imagename)}.image.pbcor'):
        logprint('Creating pbcor image {os.path.basename(imagename)}.image.pbcor')
        try:
            impbcor(imagename='{os.path.basename(imagename)}.image',
                    pbimage='{os.path.basename(imagename)}.pb',
                    outfile='{os.path.basename(imagename)}.image.pbcor',)
        except RuntimeError as ex:
            print(f"RuntimeError: {{ex}} - trying again with overwrite=False and to .image.pbcor-temp instead")
            impbcor(imagename='{os.path.basename(imagename)}.image',
                    pbimage='{os.path.basename(imagename)}.pb',
                    outfile='{os.path.basename(imagename)}.image.pbcor-temp',)
            if os.path.exists('{imagename}.image.pbcor'):
                shutil.rmtree('{imagename}.image.pbcor')
            shutil.move('{imagename}.image.pbcor-temp', '{imagename}.image.pbcor')
        print(f"Successfully created {os.path.basename(imagename)}.image.pbcor: {{os.path.exists('{imagename}.image.pbcor')}}")

    image_rbeam = ia.restoringbeam()
    image_cbeam = ia.commonbeam()
    if 'beams' in image_rbeam:
        image_manybeam = check_manybeam(image_rbeam, image_cbeam)
        assert not image_manybeam, "FAILURE: image_manybeam is True, so the conversion to convolved-beam failed"

    print("Fitsifying .image and .image.pbcor")
    exportfits(imagename='{os.path.basename(imagename)}.image',
               fitsimage='{os.path.basename(imagename)}.image.fits',
               overwrite=True
               )
    exportfits(imagename='{os.path.basename(imagename)}.image.pbcor',
               fitsimage='{os.path.basename(imagename)}.image.pbcor.fits',
               overwrite=True
               )

    print("Done with manybeam")

def check_file(fn):
    from spectral_cube import SpectralCube
    import numpy as np
    cube = SpectralCube.read(fn, format='fits' if fn.endswith('.fits') else 'casa_image', use_dask=True)
    mx = cube.max()
    if mx == 0 or np.isnan(mx):
        raise ValueError(f"File {{fn}} has a maximum value of {{mx}}")


for suffix in suffixes_to_merge_and_export:
    outfile = os.path.basename(f'{imagename}.{{suffix}}')
    if savedir and os.path.exists(savedir):
        # ensure we don't do any moving or cleanup if the files are junk
        check_file(outfile)
        check_file(outfile+".fits")

        print(f"Moving {{outfile}} to {{savedir}}")
        full_outfile = os.path.join(savedir, outfile)
        if os.path.exists(full_outfile):
            print(f"Outfile {{full_outfile}} already exists.  Check what's up.")
        elif os.path.exists(full_outfile+".fits"):
            print(f"Outfile {{full_outfile}}.fits already exists.  Check what's up.")
        else:
            shutil.move(outfile, savedir)
            shutil.move(outfile+".fits", savedir)

            if not os.path.exists(full_outfile):
                print(f"FAILURE: attempt to move {{outfile}} to {{savedir}} had no effect")
    else:
        print(f"Savedir {{savedir}} does not exist")


# Cleanup stage

print("Beginning cleanup")
for suffix in necessary_suffixes:
    if suffix not in suffixes_to_merge_and_export:
        print(f"Cleanup: Removing {imagename}.{{suffix}}", flush=True)
        shutil.rmtree(f'{imagename}.{{suffix}}')


os.chdir('{workdir}')
tclean_kwargs = {tclean_kwargs}

{rename_vis}

for vis in tclean_kwargs['vis']:
    logprint(f"Removing visibility {{vis}}")
    vis = rename_vis(vis)
    assert 'orange' not in vis
    shutil.rmtree(vis)


# if we've gotten to this point, merging has been successful so we can delete the components
print(f"Final cleanup stage: removing individual chunks", flush=True)
for suffix in necessary_suffixes:
    infiles = [f'{imagename}.{{start:04d}}.{nchan_per:03d}.{{suffix}}'
               for start in range(0, {nchan}, {nchan_per})]
    print(f"Removing {{infiles}}", flush=True)
    for fn in infiles:
        print(f"Final cleanup: Removing {{fn}}", flush=True)
        shutil.rmtree(fn)

print("Final cleanup complete", flush=True)
# report success forcefully - not clear why this is needed but CASA is not exiting cleanly
sys.exit(0)

""")

    with open(mergescriptname, 'w') as fh:
        fh.write(mergescript)


    #now = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M_%S")
    #LOGFILENAME = f"casa_log_line_{jobname}_merge_{now}.log"
    runcmd_merge = ("#!/bin/bash\n"
                    f'LOGFILENAME="{logdir}/casa_log_merge_{jobname}_${{SLURM_JOBID}}_{now}.log"\n'
                    'echo "Log file is $LOGFILENAME"\n'
                    f'xvfb-run -d /orange/adamginsburg/casa/{CASAVERSION}/bin/casa'
                    f' --nologger --nogui '
                    f' --logfile=$LOGFILENAME '
                    f' -c "execfile(\'{mergescriptname}\')"')

    slurmcmd_merge = imagename + "_slurm_cmd_merge.sh"
    with open(slurmcmd_merge, 'w') as fh:
        fh.write(runcmd_merge)

    cmd = (f'/opt/slurm/bin/sbatch --ntasks={ntasks * 4} '
           f'--mem-per-cpu={mem_per_cpu} --output={logdir}/{jobname}_merge_%j_%A_%a.log --job-name={jobname}_merge --account={account} ' +
           ('' if merge_only else f'--dependency=afterok:{jobid} ') +
           f'--qos={qos} --export=ALL --time={jobtime} {slurmcmd_merge}')

    if dry:
        print(cmd.split())
    else:
        print(cmd.split())
        sbatch = subprocess.check_output(cmd.split())
        print(f'{sbatch.decode()} with jobname={jobname}_merge')
