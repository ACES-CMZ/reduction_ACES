"""
# Serial mode (default)
python spectral_extraction_everywhere.py --serial

# Submit parallel jobs
python spectral_extraction_everywhere.py --submit-parallel --n-jobs 500

# Worker mode (called by SLURM script, not directly)
python spectral_extraction_everywhere.py --parallel-worker 0,1,2,3,4
"""
import regions
import os
import sys
import warnings
import numpy as np
from numpy import floor, ceil
from astropy import coordinates
from astropy import units as u
import glob
from spectral_cube import SpectralCube
from spectral_cube import wcs_utils
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import reproject
import scipy.ndimage

from aces import conf
basepath = conf.basepath

# Configuration
CATALOG_PATH = f'{basepath}/upload/compact_cont_source_catalog/official_cont_catalog_files/aces_compact_catalog_v0.fits'
CATALOG_NAME_PREFIX = 'ACEScatalog_v0_20260130'
PRODUCT_DIR = f"{basepath}/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/"
SPECTRUM_DIR = f"{basepath}/spectra"
FIELD_MAP_PATH = f'{basepath}/mosaics/continuum/12m_continuum_field_number_map.fits'
UID_TABLE_PATH = f'{basepath}/reduction_ACES/aces/data/tables/aces_SB_uids.csv'


class getslice(object):
    def __getitem__(self, x):
        return x


def extract_from_mask(cube, maskhdu, maskid):
    """
    """
    mask_match = maskhdu.data == maskid
    obj_slc = scipy.ndimage.find_objects(mask_match)[0]
    mask_co = mask_match[obj_slc]
    mask_ww = WCS(maskhdu.header)[obj_slc]

    mask_rep, _ = reproject.reproject_interp((mask_co, mask_ww),
                                             cube.wcs.celestial,
                                             shape_out=cube.shape[1:])

    slcs = cube.subcube_slices_from_mask(mask_rep > 0, spatial_only=True)
    scube = cube[slcs]

    spec = scube.mean(axis=(1, 2))

    return spec


def extract_background_spectrum(cube, center, source_major, source_minor,
                                inner_scale=1.5, outer_scale=3.0):
    """Extract a background annulus spectrum around a source.

    Parameters
    ----------
    cube : SpectralCube
        The spectral cube.
    center : SkyCoord
        Center of the source.
    source_major : float
        Source major axis in arcsec.
    source_minor : float
        Source minor axis in arcsec.
    inner_scale : float
        Scale factor for inner radius (relative to max(major, minor)).
    outer_scale : float
        Scale factor for outer radius (relative to max(major, minor)).

    Returns
    -------
    bg_spectrum : Spectrum1D-like
        Background spectrum (mean of annulus).
    """
    # Use circular annulus based on larger axis
    radius = max(source_major, source_minor)
    inner_radius = radius * inner_scale * u.arcsec
    outer_radius = radius * outer_scale * u.arcsec

    # Create inner and outer circular regions
    inner_reg = regions.CircleSkyRegion(center, radius=inner_radius)
    outer_reg = regions.CircleSkyRegion(center, radius=outer_radius)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # Extract outer region
        outer_scube = cube.subcube_from_regions([outer_reg])
        outer_mean = outer_scube.mean(axis=(1, 2))

        # Extract inner region
        inner_scube = cube.subcube_from_regions([inner_reg])
        inner_mean = inner_scube.mean(axis=(1, 2))

    # Calculate annulus spectrum with proper weighting
    outer_mask = outer_scube.mask.include()
    inner_mask = inner_scube.mask.include()
    n_outer = np.sum(outer_mask, axis=(1, 2)).astype(float)
    n_inner = np.sum(inner_mask, axis=(1, 2)).astype(float)
    n_annulus = n_outer - n_inner
    n_annulus[n_annulus <= 0] = 1.0

    bg_data = ((outer_mean.value * n_outer - inner_mean.value * n_inner)
               / n_annulus)

    # Create a new spectrum with the annulus data
    # We use the outer_mean as a template and create a new object with modified data
    from spectral_cube.lower_dimensional_structures import OneDSpectrum
    bg_spectrum = OneDSpectrum(value=bg_data, unit=outer_mean.unit,
                               wcs=outer_mean.wcs, meta=outer_mean.meta)

    return bg_spectrum


def load_support_data():
    """Load catalog, field map, and UID table.

    Returns
    -------
    catalog : astropy Table
        Source catalog.
    fieldmapwcs : WCS
        WCS for field map.
    fieldmapdata : array
        Field map data.
    uidtbl : astropy Table
        UID table.
    """
    catalog = Table.read(CATALOG_PATH)
    catalog.add_column([' ' * 60]*len(catalog), name='homecube')

    field_map = fits.open(FIELD_MAP_PATH)
    fieldmapwcs = WCS(field_map[0].header)
    fieldmapdata = field_map[0].data

    uidtbl = Table.read(UID_TABLE_PATH)

    return catalog, fieldmapwcs, fieldmapdata, uidtbl


def get_cube_for_source(source_row, fieldmapwcs, fieldmapdata, uidtbl):
    """Find cube files for a given source.

    Parameters
    ----------
    source_row : Table row
        Catalog row for the source.
    fieldmapwcs : WCS
        Field map WCS.
    fieldmapdata : array
        Field map data.
    uidtbl : Table
        UID table.

    Returns
    -------
    cubefns : list
        List of cube filenames.
    center : SkyCoord
        Source center coordinates.
    """
    center = SkyCoord(source_row['GLON_peak'], source_row['GLAT_peak'],
                     frame='galactic', unit=(u.deg, u.deg))

    pixcrd = list(map(lambda x: int(x), fieldmapwcs.world_to_pixel(center)))
    field_id = fieldmapdata[pixcrd[1], pixcrd[0]]

    if field_id == 0:
        return [], center

    mousid = uidtbl[field_id - 1]['12m MOUS ID']
    cubefns = glob.glob(f'{PRODUCT_DIR}/member.uid___A001_{mousid}/calibrated/working/*Sgr_A_star*.cube.I.iter1.image.pbcor.fits')

    return cubefns, center


def extract_spectra_for_source(source_row, fieldmapwcs, fieldmapdata, uidtbl,
                               spectrum_dir=SPECTRUM_DIR,
                               catalog_name_prefix=CATALOG_NAME_PREFIX,
                               skip_existing=True):
    """Extract source and background spectra for a single catalog source.

    Parameters
    ----------
    source_row : Table row
        Catalog row for the source.
    fieldmapwcs : WCS
        Field map WCS.
    fieldmapdata : array
        Field map data.
    uidtbl : Table
        UID table.
    spectrum_dir : str
        Output directory for spectra.
    catalog_name_prefix : str
        Prefix for output filenames.
    skip_existing : bool
        If True, skip sources with existing output files.

    Returns
    -------
    n_extracted : int
        Number of spectra extracted (0 or count of cubes processed).
    homecube : str
        Name of the home cube used, or empty string.
    """
    center = SkyCoord(source_row['GLON_peak'], source_row['GLAT_peak'],
                     frame='galactic', unit=(u.deg, u.deg))
    reg = regions.EllipseSkyRegion(center,
                                   width=source_row['fitted_major'] * u.arcsec,
                                   height=source_row['fitted_minor'] * u.arcsec,
                                   angle=source_row['pa'] * u.deg)

    cubefns, center = get_cube_for_source(source_row, fieldmapwcs,
                                         fieldmapdata, uidtbl)

    if not cubefns:
        return 0, ''

    n_extracted = 0
    homecube = ''

    for cubefn in cubefns:
        outfn = f"{spectrum_dir}/{catalog_name_prefix}_source{source_row['index']}_ellipseaverage_" + cubefn.split("/")[-1]

        if skip_existing and os.path.exists(outfn) and os.path.exists(outfn.replace('.fits', '_background.fits')):
            # Test that the file is openable
            fh = fits.open(outfn)
            fh.close()

            fh = fits.open(outfn.replace('.fits', '_background.fits'))
            fh.close()

            print(f"Skipping existing {outfn}", flush=True)
            homecube = os.path.basename(cubefn)
            continue

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cube = SpectralCube.read(cubefn)
            cube.allow_huge_operations = True
            ww = cube.wcs.celestial
            ww._naxis = cube.shape[1:]

        if not ww.footprint_contains(center):
            continue

        print(f"Processing source {source_row['index']} from {cubefn}", flush=True)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            scube = cube.subcube_from_regions([reg])
            avg = scube.mean(axis=(1, 2))

        hdu = avg.hdu

        y, x = scube.wcs.celestial.world_to_pixel(reg.center)
        slc = getslice()[int(floor(x)):int(ceil(x)),
                        int(floor(y)):int(ceil(y)), :]
        ww = wcs_utils.slice_wcs(scube.wcs, slc, numpy_order=False)
        hdu.header.update(ww.to_header())
        hdu.data = hdu.data[:, None, None]
        hdu.header['CATINDX'] = source_row['index']
        hdu.header['CATGLON'] = source_row['GLON_peak']
        hdu.header['CATGLAT'] = source_row['GLAT_peak']
        hdu.header['CATMAJS'] = source_row['fitted_major']
        hdu.header['CATMINS'] = source_row['fitted_minor']
        hdu.header['CATPA'] = source_row['pa']

        hdu.writeto(outfn, overwrite=True)
        print(f"  Wrote: {outfn}", flush=True)
        homecube = cubefn
        n_extracted += 1

        # Extract and save background spectrum
        bg_outfn = outfn.replace('.fits', '_background.fits')
        if not os.path.exists(bg_outfn):
            print("Extracting background spectrum...", flush=True)
            bg_spectrum = extract_background_spectrum(
                cube, center,
                source_row['fitted_major'], source_row['fitted_minor']
            )
            bg_hdu = bg_spectrum.hdu
            bg_hdu.header.update(ww.to_header())
            bg_hdu.data = bg_hdu.data[:, None, None]
            bg_hdu.header['CATINDX'] = source_row['index']
            bg_hdu.header['CATGLON'] = source_row['GLON_peak']
            bg_hdu.header['CATGLAT'] = source_row['GLAT_peak']
            bg_hdu.header['CATMAJS'] = source_row['fitted_major']
            bg_hdu.header['CATMINS'] = source_row['fitted_minor']
            bg_hdu.header['CATPA'] = source_row['pa']
            bg_hdu.header['BKGTYPE'] = 'annulus'
            bg_hdu.header['BKGINNER'] = (1.5, 'Inner radius scale factor')
            bg_hdu.header['BKGOUTER'] = (3.0, 'Outer radius scale factor')
            bg_hdu.writeto(bg_outfn, overwrite=True)
            print(f"  Background: {bg_outfn}", flush=True)

    return n_extracted, homecube


def main_serial(skip_existing=True):
    """Extract spectra for all sources serially.

    Parameters
    ----------
    skip_existing : bool
        If True, skip sources with existing output files.
    """
    print("Loading support data...", flush=True)
    catalog, fieldmapwcs, fieldmapdata, uidtbl = load_support_data()

    print(f"Processing {len(catalog)} sources...", flush=True)

    n_processed = 0
    n_with_cubes = 0

    for row in catalog:
        n_extracted, homecube = extract_spectra_for_source(
            row, fieldmapwcs, fieldmapdata, uidtbl,
            skip_existing=skip_existing
        )

        row['homecube'] = homecube

        if homecube:
            print(f"Source {row['index']} used cube {homecube}", flush=True)
            n_with_cubes += 1
        else:
            print(f"ERROR: Source {row['index']} was not in any cube", flush=True)

        n_processed += 1

        if n_processed % 100 == 0:
            print(f"Processed {n_processed}/{len(catalog)} sources, "
                  f"{n_with_cubes} with cubes", flush=True)

    n_without_cubes = len(catalog) - n_with_cubes
    print(f"\nDone! Processed {len(catalog)} sources:", flush=True)
    print(f"  {n_with_cubes} sources matched to cubes", flush=True)
    print(f"  {n_without_cubes} sources not in any cube", flush=True)


def main_parallel_worker(source_indices):
    """Worker function for parallel extraction.

    Parameters
    ----------
    source_indices : list of int
        List of source indices to process.
    """
    catalog, fieldmapwcs, fieldmapdata, uidtbl = load_support_data()

    for idx in source_indices:
        # Find the row with this index
        mask = catalog['index'] == idx
        if not mask.any():
            print(f"ERROR: Source index {idx} not found in catalog", flush=True)
            continue

        row = catalog[mask][0]
        n_extracted, homecube = extract_spectra_for_source(
            row, fieldmapwcs, fieldmapdata, uidtbl,
            skip_existing=True
        )

        if homecube:
            print(f"Source {idx} processed successfully", flush=True)
        else:
            print(f"ERROR: Source {idx} was not in any cube", flush=True)


def submit_slurm_jobs(n_jobs=500, sources_per_job=None):
    """Submit SLURM jobs for parallel extraction.

    Parameters
    ----------
    n_jobs : int
        Number of SLURM jobs to create.
    sources_per_job : int, optional
        Number of sources per job. If None, divides sources evenly.
    """
    catalog = Table.read(CATALOG_PATH)
    n_sources = len(catalog)
    source_indices = catalog['index'].data

    if sources_per_job is None:
        sources_per_job = int(np.ceil(n_sources / n_jobs))

    print(f"Submitting {n_jobs} jobs for {n_sources} sources "
          f"({sources_per_job} sources/job)", flush=True)

    # Create job script template
    # 4h is not enough for 100 spectra
    # 24h is hopefully overkill....
    job_script_template = """#!/bin/bash
#SBATCH --job-name=aces_spec_extract_{job_id}
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_spec_extract_%j_{job_id}.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH --qos=astronomy-dept-b
#SBATCH --account=astronomy-dept

source /orange/adamginsburg/miniconda3/bin/activate /blue/adamginsburg/adamginsburg/miniconda3/envs/python313

cd /blue/adamginsburg/adamginsburg/ACES/logs

python /orange/adamginsburg/ACES/reduction_ACES/aces/analysis/spectral_extraction_everywhere.py --parallel-worker {source_list}
"""

    # Create log directory
    log_dir = f"/blue/adamginsburg/adamginsburg/ACES/logs"
    os.makedirs(log_dir, exist_ok=True)

    # Split sources into chunks
    for job_id in range(n_jobs):
        start_idx = job_id * sources_per_job
        end_idx = min((job_id + 1) * sources_per_job, n_sources)

        if start_idx >= n_sources:
            break

        job_sources = source_indices[start_idx:end_idx]
        source_list = ','.join(map(str, job_sources))

        # Create job script
        job_script = job_script_template.format(
            job_id=job_id,
            basepath=basepath,
            source_list=source_list
        )

        script_path = f"{log_dir}/spec_extract_job_{job_id}.sh"
        with open(script_path, 'w') as f:
            f.write(job_script)

        # Submit job
        os.system(f"sbatch {script_path}")
        print(f"Submitted job {job_id} for sources {start_idx}-{end_idx}", flush=True)

    print(f"\nSubmitted {job_id + 1} jobs", flush=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Extract spectra from ACES cubes for catalog sources'
    )
    parser.add_argument('--serial', action='store_true',
                       help='Run in serial mode (default)')
    parser.add_argument('--submit-parallel', action='store_true',
                       help='Submit SLURM jobs for parallel processing')
    parser.add_argument('--parallel-worker', type=str,
                       help='Run as parallel worker with comma-separated source indices')
    parser.add_argument('--n-jobs', type=int, default=500,
                       help='Number of SLURM jobs (for --submit-parallel)')
    parser.add_argument('--no-skip-existing', action='store_true',
                       help='Re-extract even if output files exist')

    args = parser.parse_args()

    if args.parallel_worker:
        # Parse source indices
        source_indices = [int(x) for x in args.parallel_worker.split(',')]
        print(f"Running parallel worker for {len(source_indices)} sources", flush=True)
        main_parallel_worker(source_indices)
    elif args.submit_parallel:
        submit_slurm_jobs(n_jobs=args.n_jobs)
    else:
        # Default: serial mode
        main_serial(skip_existing=not args.no_skip_existing)
