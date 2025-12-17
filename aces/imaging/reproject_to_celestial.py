"""
Reproject final products to RA/Dec.

Target header is header_12m_radec.hdr

target folder is {ACES_rootdir}/products_for_ALMA/group.uid___A001_X1590_X30a9.lp_slongmore/.  All non-PV files in this directory should be reprojected.

renaming: reprojected files should go in {ACES_rootdir}//products_for_ALMA/group.uid___A001_X1590_X30a9.lp_slongmore/icrs/, and .icrs should be added to each file's name after `cmz_mosaic`

sbatch --job-name=reproject_aces_cel --account=astronomy-dept --qos=astronomy-dept-b --nodes=1 --ntasks=32 --mem=512gb --time=96:00:00 --output=/blue/adamginsburg/adamginsburg/logs/aces_reproject_celestial_%j.log --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python312/bin/python /orange/adamginsburg/ACES/reduction_ACES/aces/imaging/reproject_to_celestial.py --serial"


To use dask, we need to handle:
https://github.com/dask/dask/issues/11850#issuecomment-3035792958
https://github.com/dask/dask/pull/11524
"""

import os
import tempfile
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
import numpy as np
import warnings
from tqdm import tqdm
import dask.array as da
from dask.diagnostics import ProgressBar
import dask
from dask import delayed

# Configuration
ACES_rootdir = Path("/orange/adamginsburg/ACES")
header_file = ACES_rootdir / "reduction_ACES/aces/imaging/data/header_12m_radec.hdr"
source_dir = ACES_rootdir / "products_for_ALMA/group.uid___A001_X1590_X30a9.lp_slongmore"
output_dir = source_dir / "icrs"

# Set up Dask temporary directory for Slurm environment
# Use SLURM_TMPDIR if available, otherwise fall back to /blue scratch space
if 'SLURM_TMPDIR' in os.environ:
    dask_tempdir = Path(os.environ['SLURM_TMPDIR']) / 'dask-worker-space'
else:
    dask_tempdir = Path('/blue/adamginsburg/adamginsburg/tmp/dask-worker-space')

dask_tempdir.mkdir(parents=True, exist_ok=True)

# Configure Dask
dask.config.set(temporary_directory=str(dask_tempdir))
dask.config.set({'array.slicing.split_large_chunks': True})
print(f"Dask temporary directory: {dask_tempdir}")


def load_target_header(header_path):
    """Load the target header from a .hdr file."""
    print(f"Loading target header from {header_path}")
    header = fits.Header.fromtextfile(str(header_path))
    return header


def get_files_to_reproject(directory):
    """
    Get all FITS files in the directory, excluding PV files and files without 'cmz_mosaic'.

    Returns:
        List of Path objects for files to reproject
    """
    directory = Path(directory)
    all_fits = list(directory.glob("*.fits"))

    # Exclude PV files and files without cmz_mosaic in their name
    valid_files = [f for f in all_fits if ".PV_" not in f.name and "cmz_mosaic" in f.name]

    print(f"Found {len(all_fits)} total FITS files")
    print(f"Excluding {len(all_fits) - len(valid_files)} files (PV files or missing 'cmz_mosaic')")
    print(f"Will reproject {len(valid_files)} files")

    return sorted(valid_files)


@delayed
def load_fits_slice(filepath, slices):
    """Lazily load a slice from a FITS file using memmap."""
    with fits.open(filepath, memmap=True, mode='readonly') as hdul:
        return np.array(hdul[0].data[slices], dtype=np.float32)


def reproject_file(input_path, target_header, output_path, block_size='auto'):
    """
    Reproject a FITS file to the target coordinate system using Dask arrays and memory-mapped output.

    Parameters:
        input_path: Path to input FITS file
        target_header: Target header for reprojection
        output_path: Path for output FITS file
        block_size: Chunk size for Dask arrays (spectral axis for cubes)
    """
    print(f"\nReprojecting: {input_path.name}")

    # Read header only (no data loading)
    with fits.open(input_path, memmap=True, mode='readonly') as hdul:
        hdu = hdul[0]
        input_header = hdu.header.copy()

        # Prepare header for output
        new_header = input_header.copy()
        new_header.update(target_header)

        # Get shape from header (not data) to avoid loading
        naxis = input_header.get('NAXIS', 0)
        if naxis == 2:
            input_shape = (input_header['NAXIS2'], input_header['NAXIS1'])
            output_shape = (int(target_header['NAXIS2']), int(target_header['NAXIS1']))
        elif naxis == 3:
            input_shape = (input_header['NAXIS3'], input_header['NAXIS2'], input_header['NAXIS1'])
            output_shape = (input_header['NAXIS3'], int(target_header['NAXIS2']), int(target_header['NAXIS1']))
            if block_size == 'auto':
                chunk_size = min(50, input_shape[0])
            else:
                chunk_size = block_size
            chunks = (chunk_size, input_shape[1], input_shape[2])
        elif naxis == 4:
            input_shape = (input_header['NAXIS4'], input_header['NAXIS3'], input_header['NAXIS2'], input_header['NAXIS1'])
            output_shape = (input_header['NAXIS4'], input_header['NAXIS3'], int(target_header['NAXIS2']), int(target_header['NAXIS1']))
            if block_size == 'auto':
                chunk_size = min(50, input_shape[1])
            else:
                chunk_size = block_size
            chunks = (1, chunk_size, input_shape[2], input_shape[3])
        else:
            raise ValueError(f"Unsupported NAXIS={naxis}")

        # Create lazy Dask array using delayed loading
        if naxis > 2:
            print(f"    Creating lazy Dask array with shape {input_shape}, chunks {chunks}")
            # Build Dask array from delayed FITS file loads
            lazy_arrays = []
            if naxis == 3:
                for i in range(input_shape[0]):
                    lazy_slice = da.from_delayed(
                        load_fits_slice(str(input_path), (i, slice(None), slice(None))),
                        shape=(input_shape[1], input_shape[2]),
                        dtype=np.float32
                    )
                    lazy_arrays.append(lazy_slice)
                input_dask = da.stack(lazy_arrays, axis=0)
            else:
                raise NotImplementedError("4D lazy loading not yet implemented")
        else:
            # For 2D, create a lazy delayed load with float32 conversion
            lazy_2d = da.from_delayed(
                load_fits_slice(str(input_path), (slice(None), slice(None))),
                shape=(input_shape[0], input_shape[1]),
                dtype=np.float32
            )
            input_dask = lazy_2d

    # Create output FITS file with memory-mapped array
    print(f"    Creating memory-mapped output with shape {output_shape}")
    output_hdu = fits.PrimaryHDU(data=np.zeros(output_shape, dtype=np.float32), header=new_header)
    output_hdu.writeto(output_path, overwrite=True)

    # Re-open with memmap for writing
    with fits.open(output_path, mode='update', memmap=True) as output_hdul:
        # Reproject into the memory-mapped array
        if naxis == 2:
            print(f"    Reprojecting 2D image")
            reproject_interp((input_dask, WCS(input_header)), target_header,
                            output_array=output_hdul[0].data, return_footprint=False,
                            parallel=True)
        elif naxis == 3:
            print(f"    Reprojecting 3D cube with Dask")
            n_planes = output_shape[0]
            wcs_3d = WCS(input_header)
            wcs_2d = wcs_3d.celestial

            for i in tqdm(range(n_planes), desc="    Planes", leave=False):
                # Dask slice is lazy - only computed when needed
                plane_data = input_dask[i, :, :]
                reproject_interp((plane_data, wcs_2d), target_header,
                                output_array=output_hdul[0].data[i, :, :],
                                return_footprint=False, parallel=True)
        else:
            # really this is a 'NotRelevantError'
            raise NotImplementedError("Reprojection for 4D data is not implemented yet.")

        # Flush changes to disk
        output_hdul.flush()

    print(f"    ✓ Complete")


def calculate_memory_allocation(file_size_bytes):
    """
    Calculate memory allocation based on file size.
    Allocate 2x file size, rounded up to nearest power-of-2-like value.

    Parameters:
        file_size_bytes: File size in bytes

    Returns:
        Memory allocation in GB (256, 384, or 512)
    """
    # Convert to GB and multiply by 2
    file_size_gb = file_size_bytes / (1024**3)
    needed_gb = file_size_gb * 2

    # Round up to standard memory allocations
    if needed_gb <= 256:
        return 256
    elif needed_gb <= 384:
        return 384
    else:
        return 512


def main_serial():
    """Serial reprojection workflow (original main function)."""
    print("ACES Data Reprojection to ICRS (RA/Dec) - Serial Mode")

    # Load target header
    target_header = load_target_header(header_file)

    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True, parents=True)
    print(f"\nOutput directory: {output_dir}")

    # Get files to reproject
    files_to_reproject = get_files_to_reproject(source_dir)

    # Reproject each file
    print(f"Starting reprojection of {len(files_to_reproject)} files")

    for input_path in files_to_reproject:
        # Insert .icrs after cmz_mosaic in the filename
        output_path = output_dir / input_path.name.replace('cmz_mosaic.', 'cmz_mosaic.icrs.')

        if os.path.exists(output_path):
            print(f"    Skipping {output_path.name}, already exists.")
        else:
            #reproject_file(input_path, target_header, output_path)
            reproject_interp(input_path, target_header, output_path)

    # Summary
    print(f"Reprojection complete!  Total: {len(files_to_reproject)}")


def main():
    """Submit SLURM array job for parallel reprojection."""
    import subprocess
    import sys

    print("ACES Data Reprojection to ICRS (RA/Dec) - SLURM Array Mode")

    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True, parents=True)
    print(f"Output directory: {output_dir}")

    # Get files to reproject
    files_to_reproject = get_files_to_reproject(source_dir)

    # Filter out already completed files
    files_needed = []
    for input_path in files_to_reproject:
        output_path = output_dir / input_path.name.replace('cmz_mosaic.', 'cmz_mosaic.icrs.')
        if not os.path.exists(output_path):
            files_needed.append(input_path)
        else:
            print(f"    Skipping {output_path.name}, already exists.")

    if not files_needed:
        print("\nNo files need reprojection!")
        return

    print(f"\nWill submit {len(files_needed)} jobs to SLURM array")

    # Create a file list for the array job to read
    filelist_path = output_dir / "reproject_filelist.txt"
    with open(filelist_path, 'w') as f:
        for fpath in files_needed:
            f.write(f"{fpath}\n")

    print(f"File list written to: {filelist_path}")

    # Calculate memory allocations for each file and write to file
    memory_allocations = []
    memory_file_path = output_dir / "reproject_memory_allocations.txt"
    with open(memory_file_path, 'w') as f:
        for input_path in files_needed:
            file_size = input_path.stat().st_size
            mem_gb = calculate_memory_allocation(file_size)
            memory_allocations.append(mem_gb)
            f.write(f"{mem_gb}\n")
            print(f"  {input_path.name}: {file_size/(1024**3):.2f} GB -> allocate {mem_gb} GB")

    print(f"Memory allocations written to: {memory_file_path}")

    # Submit individual SLURM jobs with appropriate memory for each file
    script_path = Path(__file__).absolute()
    python_exe = sys.executable

    job_ids = []
    for i, (input_path, mem_gb) in enumerate(zip(files_needed, memory_allocations)):
        sbatch_command = [
            'sbatch',
            f'--job-name=reproj_cel_{i}',
            '--account=astronomy-dept',
            '--qos=astronomy-dept-b',
            '--nodes=1',
            '--ntasks=16',
            f'--mem={mem_gb}gb',
            '--time=96:00:00',
            f'--output=/blue/adamginsburg/adamginsburg/logs/aces_reproject_celestial_{i}_%j.log',
            '--wrap',
            f'{python_exe} {script_path} --array-task --task-id={i}'
        ]

        result = subprocess.run(sbatch_command, capture_output=True, text=True)

        if result.returncode == 0:
            job_id = result.stdout.strip().split()[-1]
            job_ids.append(job_id)
            print(f"  Submitted job {i}/{len(files_needed)}: {job_id} ({mem_gb}GB) - {input_path.name}")
        else:
            print(f"  ✗ Error submitting job {i}: {result.stderr}")

    print(f"\n✓ Submitted {len(job_ids)} individual SLURM jobs")

    return 0


def array_task_worker():
    """Worker function for SLURM array job - processes one file."""
    import sys

    # Get task ID from command line or environment
    task_id = None
    for arg in sys.argv:
        if arg.startswith('--task-id='):
            task_id = int(arg.split('=')[1])
            break

    if task_id is None and 'SLURM_ARRAY_TASK_ID' in os.environ:
        task_id = int(os.environ['SLURM_ARRAY_TASK_ID'])

    if task_id is None:
        print("Error: No task ID found in arguments or environment")
        sys.exit(1)

    print(f"Task {task_id} starting")

    # Read file list
    filelist_path = output_dir / "reproject_filelist.txt"
    with open(filelist_path, 'r') as f:
        files = [Path(line.strip()) for line in f]

    if task_id >= len(files):
        print(f"Error: Task ID {task_id} is out of range (only {len(files)} files)")
        sys.exit(1)

    input_path = files[task_id]
    output_path = output_dir / input_path.name.replace('cmz_mosaic.', 'cmz_mosaic.icrs.')

    print(f"Processing file {task_id + 1}/{len(files)}: {input_path.name}")

    # Load target header
    target_header = load_target_header(header_file)

    # Reproject this file
    try:
        reproject_file(input_path, target_header, output_path)
        print(f"✓ Task {task_id} completed successfully")
    except Exception as e:
        print(f"✗ Task {task_id} failed with error:")
        print(e)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    import sys

    # Check if this is an array task worker or the main submission script
    if '--array-task' in sys.argv or 'SLURM_ARRAY_TASK_ID' in os.environ:
        array_task_worker()
    elif '--serial' in sys.argv:
        main_serial()
    else:
        main()