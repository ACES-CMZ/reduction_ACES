"""
Processes files in "products_for_ALMA" directory to ensure compliance with ALMA Science Archive requirements.
- Standardises mandatory FITS keywords (OBJECT, TELESCOP, INSTRUME, RADESYS, ORIGIN)
- Converts velocity axes to frequency axes where needed
- Adds degenerate frequency and/or Stokes axes where needed
- Transfers frequency WCS information from parent cubes to moment maps
- Standardises BUNIT to 'Jy/beam'
- Computes and populates DATAMIN/DATAMAX keywords
- Adds a new TARGET keyword for non-mosaic products, which is the central Galactic coordinate of the field
- Adds a standard DATE-OBS keyword for all mosaic products (default to ACES start date)
- Adds a standard ORIGIN keyword (links to ACES GitHub)
"""

import re
import os
import sys
import time
import warnings
from pathlib import Path
from typing import Optional
import numpy as np
import dask
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import add_stokes_axis_to_wcs
from spectral_cube import SpectralCube
from dask.diagnostics import ProgressBar

BASE_DIR = Path("/orange/adamginsburg/ACES/products_for_ALMA/")
STATS_DIR = Path("/orange/adamginsburg/ACES/reduction_ACES/aces/data/tables/")
CUBE_STATS_PATH = STATS_DIR / "cube_stats.ecsv"
FEATHERED_STATS_PATH = STATS_DIR / "feathered_cube_stats.ecsv"

DEFAULT_DATE_OBS = "2021-10-17T01:21:04"
DEFAULT_ORIGIN = "ACES team: see https://github.com/ACES-CMZ/"

# Global cache for stats tables
_CUBE_STATS = None
_FEATHERED_STATS = None

GAL_COORD_RE = re.compile(r"(G\d{1,3}\.\d+[pm]\d{1,2}\.\d+)")

CONTINUUM_PATTERNS = {
    ".cont.93.7GHz_bw4.7GHz": {"freq": 93.7e9, "bw": 15.5e9},
    ".cont.99.6GHz_bw3.8GHz": {"freq": 99.6e9, "bw": 3.8e9},
    ".cont.86.6GHz_bw1.0GHz": {"freq": 86.6e9, "bw": 1.2e9},
}

MOMENT_SUFFIXES = (
    ".integrated_intensity.fits", ".mad_std.fits",
    ".peak_intensity.fits", ".velocity_at_peak_intensity.fits",
)

MALFORMED_BUNIT = {'Jy beam-1', 'Jy.beam-1', 'beam-1.Jy', 'Jy beam^-1', 'Jy.beam^-1'}


def has_datamax(fits_path: Path) -> bool:
    """Check if a FITS file has DATAMAX header using astropy.io.fits."""
    try:
        with fits.open(fits_path, mode='readonly', memmap=True) as hdul:
            return 'DATAMAX' in hdul[0].header
    except (OSError, IOError):
        # If file can't be opened, there's a big problem
        raise


def load_stats_tables():
    """Load the stats tables once and cache them."""
    global _CUBE_STATS, _FEATHERED_STATS

    if _CUBE_STATS is None and CUBE_STATS_PATH.exists():
        print(f"Loading cube stats from {CUBE_STATS_PATH}", flush=True)
        _CUBE_STATS = Table.read(CUBE_STATS_PATH)

    if _FEATHERED_STATS is None and FEATHERED_STATS_PATH.exists():
        print(f"Loading feathered stats from {FEATHERED_STATS_PATH}", flush=True)
        _FEATHERED_STATS = Table.read(FEATHERED_STATS_PATH)

    return _CUBE_STATS, _FEATHERED_STATS


def get_original_filename(current_path: Path) -> str:
    """
    Map from current filename to original filename used in stats tables.

    The files were renamed by update_fns.py which changed:
    'lp_2021.1.00172.L.slongmore' -> 'lp_slongmore'

    Also need to handle the different naming conventions:
    - products_for_ALMA uses: group.uid___A001_X1590_X30a9.lp_slongmore.*.fits
    - cube_stats uses: uid___*.image
    - feathered_stats uses: Sgr_A_st*.fits
    """
    filename = current_path.name

    # Try to extract the core filename for matching
    # For feathered files: look for pattern like Sgr_A_st_*.TP_7M_12M_feather*
    if '12m7mTP' in filename or 'feather' in filename.lower():
        # Try to match feathered naming pattern
        # e.g., group...lp_slongmore.cmz_mosaic.icrs.12m7mTP.CH3CHO.cube.pbcor.fits
        # might match: Sgr_A_st_*.TP_7M_12M_feather_all.SPW_*.image.statcont.contsub.fits
        # This is complex, so we'll try to extract key parts
        return filename  # For now return as-is, will refine
    else:
        # For non-feathered: might be in cube_stats
        # Look for uid pattern in the filename
        # e.g., group.uid___A001_X1590_X30a9.lp_slongmore.G000.029p0.041.12m7m.89.1-89.2GHz.cube.pbcor.fits
        # should map to something in cube_stats
        return filename


def lookup_precomputed_stats(fits_path: Path) -> Optional[tuple[float, float]]:
    """
    Look up precomputed datamin and datamax from stats tables.

    NOTE: Currently disabled - the mapping between current filenames in products_for_ALMA
    and the original filenames in the stats tables is complex:
    1. Stats tables (cube_stats, feathered_cube_stats) were created from intermediate
       processing files with different UIDs and naming conventions
    2. Files in products_for_ALMA went through renaming (update_fns.py) and reorganization
    3. No direct mapping file exists between the two naming schemes

    TODO: Create a proper mapping from products_for_ALMA filenames back to the original
    intermediate filenames used when generating the stats tables. This would require:
    - Tracking the processing pipeline to link final products to intermediate cubes
    - Building a lookup table during the processing/renaming steps
    - Or recomputing stats from the final products and saving a new stats table

    For now, this function returns None and datamin/datamax are computed fresh.

    Returns
    -------
    tuple of (datamin, datamax) or None if not found (currently always None)
    """
    # Disabled for now - uncomment and fix matching logic when mapping is available
    return None

    # The code below shows the intended structure:
    #
    # cube_stats, feathered_stats = load_stats_tables()
    #
    # if cube_stats is None and feathered_stats is None:
    #     return None
    #
    # filename = fits_path.name
    #
    # # Extract UID and SPW for matching
    # uid_match = re.search(r'(uid___A\d+_X[0-9a-f]+_X[0-9a-f]+)', filename)
    # uid_str = uid_match.group(1) if uid_match else None
    #
    # if cube_stats is not None and uid_str:
    #     for row in cube_stats:
    #         table_filename = row['filename']
    #         if uid_str in table_filename:
    #             spw_match = re.search(r'spw(\d+)', filename.lower())
    #             table_spw_match = re.search(r'spw(\d+)', table_filename.lower())
    #
    #             if spw_match and table_spw_match and spw_match.group(1) == table_spw_match.group(1):
    #                 datamin = float(row['min'].value if hasattr(row['min'], 'value') else row['min'])
    #                 datamax = float(row['max'].value if hasattr(row['max'], 'value') else row['max'])
    #                 return (datamin, datamax)
    #
    # return None


def extract_gal_target(name: str) -> Optional[str]:
    return m.group(1) if (m := GAL_COORD_RE.search(name)) else None


def get_axis_types(header) -> Optional[list]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return WCS(header, naxis=header.get('NAXIS')).world_axis_physical_types


def axis_is_spatial(t: str) -> bool:
    return bool(t) and (t.startswith('pos.') or t.startswith('celestial.'))


def axis_is_spectral(t: str) -> bool:
    return bool(t) and (t.startswith('em.') or t.startswith('spect.') or
                        'freq' in t.lower() or 'velo' in t.lower())


def is_velocity_axis(header) -> tuple[bool, Optional[str]]:
    """Check if axis 3 is velocity."""
    if header.get("NAXIS", 0) < 3:
        return False, None
    ctype3 = header.get("CTYPE3", "").upper()
    if "VELO" not in ctype3 and "VRAD" not in ctype3:
        return False, None
    if not header.get("RESTFRQ"):
        return False, "Velocity axis but missing RESTFRQ"
    return True, None


def convert_velocity_to_frequency(header) -> list[str]:
    """Convert axis 3 from velocity to frequency."""
    restfrq = header.get("RESTFRQ")
    cunit = header.get("CUNIT3", "m/s")
    rest_freq = restfrq * u.Hz
    vel_ref = header["CRVAL3"] * u.Unit(cunit)
    vel_delta = header["CDELT3"] * u.Unit(cunit)
    equiv = u.doppler_radio(rest_freq)
    freq_ref = vel_ref.to(u.Hz, equivalencies=equiv)
    freq_delta = (vel_ref + vel_delta).to(u.Hz, equivalencies=equiv) - freq_ref

    header["CTYPE3"] = "FREQ"
    header["CUNIT3"] = "Hz"
    header["CRVAL3"] = freq_ref.value
    header["CDELT3"] = freq_delta.value
    header["SPECSYS"] = "LSRK"
    header["VELREF"] = 257

    return ["converted axis 3 velocity -> frequency; set SPECSYS='LSRK', VELREF=257"]


def find_cube_for_moment(path: Path) -> Optional[Path]:
    """Find corresponding cube file for a moment map."""
    name_lower = path.name.lower()
    for suffix in MOMENT_SUFFIXES:
        if name_lower.endswith(suffix.lower()):
            prefix = path.name[:len(path.name) - len(suffix)]
            for cube_suffix in (".cube.pbcor.fits", ".cube.downsampled_spectrally.pbcor.fits",
                                ".cube.downsampled_spatially.pbcor.fits"):
                if (cube := path.parent / f"{prefix}{cube_suffix}").exists():
                    return cube


def extract_spectral_wcs(cube_path: Path) -> Optional[dict]:
    """Extract spectral axis WCS info from a cube."""
    with fits.open(cube_path, mode='readonly') as hdul:
        hdr = hdul[0].header
        naxis = hdr.get("NAXIS")
        if not isinstance(naxis, int) or naxis < 3:
            raise ValueError(f"Cube found at {cube_path} but NAXIS is {naxis} (expected >= 3)")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wcs = WCS(hdr, naxis=naxis)
            axis_types = wcs.world_axis_physical_types

            spec_idx = next((i for i, t in enumerate(axis_types) if axis_is_spectral(t)), None)
            if spec_idx is None:
                raise ValueError(f"Cube found at {cube_path} but no spectral axis identified")

            ax = spec_idx + 1
            return {
                "CTYPE3": hdr.get(f"CTYPE{ax}"),
                "CRVAL3": hdr.get(f"CRVAL{ax}"),
                "CDELT3": hdr.get(f"CDELT{ax}"),
                "CRPIX3": hdr.get(f"CRPIX{ax}"),
                "CUNIT3": hdr.get(f"CUNIT{ax}", "Hz"),
                "RESTFRQ": hdr.get("RESTFRQ"),
            }


def identify_continuum(path: Path) -> Optional[dict]:
    return next((p for pat, p in CONTINUUM_PATTERNS.items() if pat in path.name), None)


def should_add_freq_axis(header, cont_params: Optional[dict]) -> bool:
    if not cont_params or header.get("NAXIS") != 2:
        return False
    axis_types = get_axis_types(header)
    return (axis_types and len(axis_types) >= 2 and
            axis_is_spatial(axis_types[0]) and axis_is_spatial(axis_types[1]))


def should_add_stokes_axis(header, is_pv: bool) -> bool:
    if is_pv or header.get("NAXIS") != 3:
        return False
    axis_types = get_axis_types(header)
    return (axis_types and len(axis_types) >= 3 and
            axis_is_spatial(axis_types[0]) and axis_is_spatial(axis_types[1]) and
            axis_is_spectral(axis_types[2]))


def apply_fixes(fits_path: Path) -> tuple[bool, list[str]]:
    """Apply all fixes to a FITS file. Returns (changed, error_messages)."""
    errors = []
    fname_lower = fits_path.name.lower()
    t0 = time.time()

    with fits.open(fits_path, mode="update") as hdul:

        print(f"Loaded header of {os.path.basename(fits_path)}. dt={time.time() - t0:.2f}s", flush=True)
        h0 = hdul[0]
        hdr = h0.header

        # Identify file type
        is_model_resid = ".model" in fname_lower or ".residual" in fname_lower
        is_pv = ".pv" in fname_lower
        is_mosaic = "cmz_mosaic" in fname_lower
        is_moment = any(fname_lower.endswith(s.lower()) for s in MOMENT_SUFFIXES)

        target = None if is_mosaic else extract_gal_target(fits_path.name)
        cont_params = identify_continuum(fits_path)

        # Handle moment maps
        cube_wcs = None
        if is_moment and (cube_file := find_cube_for_moment(fits_path)):
            print("Extracting spectral WCS from cube:", cube_file, flush=True)
            cube_wcs = extract_spectral_wcs(cube_file)

        # === Apply fixes ===

        # 1. Static mandatory keywords
        hdr["OBJECT"] = "Sgr_A_star"
        hdr["TELESCOP"] = "ALMA"
        hdr["INSTRUME"] = "ALMA"
        hdr["RADESYS"] = "ICRS"

        # ORIGIN / CASAVERS
        cur_origin = hdr.get("ORIGIN")
        if cur_origin and "casa" in str(cur_origin).lower():
            hdr["CASAVERS"] = (cur_origin, "CASA version used for processing")
        hdr["ORIGIN"] = DEFAULT_ORIGIN

        # DATE-OBS
        if not (dob := hdr.get("DATE-OBS")):
            hdr["DATE-OBS"] = DEFAULT_DATE_OBS
        elif " " in str(dob) and "T" not in str(dob):
            hdr["DATE-OBS"] = str(dob).replace(" ", "T", 1)

        if target:
            hdr["TARGET"] = target

        print(f"Loading cube from {fits_path} for further analysis... dt={time.time() - t0:.2f}s", flush=True)
        # 2. Dimension updates
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                cube = SpectralCube.read(fits_path, use_dask=True)
        except Exception as ex:
            print("This should only happen for 2D images and so we shouldn't be catching generic exceptions")
            print(ex, flush=True)
            cube = h0.data

        data = h0.data

        # Moment maps: 2D -> 3D
        if is_moment and cube_wcs and hdr.get("NAXIS") == 2:
            if should_add_freq_axis(hdr, {"freq": 1, "bw": 1}) and data is not None:
                print(f"Converting moment map from 2d to 3d. dt={time.time() - t0:.2f}s", flush=True)
                h0.data = data[np.newaxis, ...]
                data = h0.data
                hdr["NAXIS"] = 3
                hdr["NAXIS3"] = 1
                hdr["CTYPE3"] = cube_wcs["CTYPE3"]
                hdr["CRVAL3"] = float(cube_wcs["CRVAL3"])
                hdr["CRPIX3"] = float(cube_wcs["CRPIX3"])
                hdr["CDELT3"] = float(cube_wcs["CDELT3"])
                hdr["CUNIT3"] = cube_wcs["CUNIT3"]
                if cube_wcs.get("RESTFRQ"):
                    hdr["RESTFRQ"] = cube_wcs["RESTFRQ"]

        # Continuum: 2D -> 3D
        elif should_add_freq_axis(hdr, cont_params) and data is not None:
            print(f"Converting continuum map from 2d to 3d. dt={time.time() - t0:.2f}s", flush=True)
            h0.data = data[np.newaxis, ...]
            data = h0.data
            hdr["NAXIS"] = 3
            hdr["NAXIS3"] = 1
            hdr["CTYPE3"] = "FREQ"
            hdr["CRVAL3"] = float(cont_params['freq'])
            hdr["CRPIX3"] = 1.0
            hdr["CDELT3"] = float(cont_params['bw'])
            hdr["CUNIT3"] = "Hz"

        # 3D -> 4D with Stokes
        if should_add_stokes_axis(hdr, is_pv) and data is not None:
            print(f"Adding Stokes axis to make 4D. dt={time.time() - t0:.2f}s", flush=True)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                old_wcs = WCS(hdr, naxis=3)
                new_wcs = add_stokes_axis_to_wcs(old_wcs, old_wcs.pixel_n_dim)
                h0.data = data[np.newaxis, ...]
                data = h0.data
                hdr.update(new_wcs.to_header())
                hdr["NAXIS"] = 4
                hdr["NAXIS4"] = 1
                hdr["CUNIT4"] = ""

        # Ensure CUNIT4 on existing 4D files
        if hdr.get("NAXIS") == 4 and "CUNIT4" not in hdr:
            hdr["CUNIT4"] = ""

        # 3. Velocity -> Frequency conversion
        if is_velocity_axis(hdr)[0]:
            convert_velocity_to_frequency(hdr)
        else:
            hdr["SPECSYS"] = "LSRK"
            hdr["VELREF"] = 257

        # 4. Data min/max
        if data is not None and ('DATAMIN' not in hdr or 'DATAMAX' not in hdr):
            # First try to get precomputed values from stats tables
            precomputed = lookup_precomputed_stats(fits_path)

            if precomputed is not None:
                datamin, datamax = precomputed
                print(f"Using precomputed datamin/datamax from stats tables. dt={time.time() - t0:.2f}s", flush=True)
            else:
                # Calculate from scratch
                print(f"Calculating datamin/datamax for data with size {data.size/1024**3}Gpix... dt={time.time() - t0:.2f}s", flush=True)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    try:
                        with dask.config.set(scheduler='threads'):
                            with ProgressBar():
                                datamax = cube.max()
                                datamin = cube.min()

                        # this shouldn't be necessary, but it has happened: they should _always_ be qty's
                        if hasattr(datamax, 'value'):
                            datamax = datamax.value
                        if hasattr(datamin, 'value'):
                            datamin = datamin.value
                    except TypeError:
                        # it's a moment map, we can afford to load the whole data
                        datamax = np.nanmax(cube)
                        datamin = np.nanmin(cube)

            if np.isnan(datamax) or np.isnan(datamin):
                errors.append("Could not compute DATAMIN/DATAMAX (resulted in NaN)")
                print(f"Warning: could not compute DATAMIN/DATAMAX for file {fits_path} (resulted in NaN). dt={time.time() - t0:.2f}s", flush=True)
            else:
                if "DATAMIN" not in hdr:
                    hdr["DATAMIN"] = float(datamin)
                if "DATAMAX" not in hdr:
                    hdr["DATAMAX"] = float(datamax)
            print(f"Computed datamin/datamax for file {fits_path}, datamax={hdr.get('DATAMAX')}, datamin={hdr.get('DATAMIN')}. dt={time.time() - t0:.2f}s", flush=True)
        else:
            print(f"for file {fits_path}, datamax={hdr.get('DATAMAX')}, datamin={hdr.get('DATAMIN')}. dt={time.time() - t0:.2f}s", flush=True)

        # 5. BUNIT standardisation
        if not is_model_resid:
            bunit = hdr.get("BUNIT")
            if bunit and bunit.strip() in MALFORMED_BUNIT:
                hdr["BUNIT"] = "Jy/beam"

        print(f"Starting final flush to disk... dt={time.time() - t0:.2f}s", flush=True)
        hdul.flush()
        return True, errors


def main(force: bool = False):
    """
    Process FITS files to update headers.

    Parameters
    ----------
    force : bool, optional
        If True, process all files even if they already have DATAMAX.
        If False (default), only process files missing DATAMAX.
    """
    fits_list = sorted([
        f for f in BASE_DIR.glob("*/*.fits")
    ])
    # PV files aren't valid anyway, no need to update them
    fits_list = [x for x in fits_list if 'PV' not in str(x)]

    if not force:
        print(f"Found {len(fits_list)} FITS files (excluding PV), checking for DATAMAX...", flush=True)
        # Filter out files that already have DATAMAX
        fits_list = [f for f in fits_list if not has_datamax(f)]
        print(f"Total of {len(fits_list)} FITS files without DATAMAX in {BASE_DIR}", flush=True)
    else:
        print(f"Force mode: processing all {len(fits_list)} FITS files in {BASE_DIR}", flush=True)

    total = len(fits_list)

    # Check if running as SLURM array job
    if os.getenv('SLURM_ARRAY_TASK_ID') is not None:
        slurm_array_task_id = int(os.getenv('SLURM_ARRAY_TASK_ID'))

        if slurm_array_task_id >= total:
            print(f"SLURM_ARRAY_TASK_ID={slurm_array_task_id} >= total files ({total}), nothing to do", flush=True)
            return

        # Process only the file assigned to this array task
        f = fits_list[slurm_array_task_id]
        print(f"SLURM array task {slurm_array_task_id}: Processing {f}", flush=True)
        ok, errs = apply_fixes(f)

        if ok:
            print(f"Successfully updated: {f.name}", flush=True)
        else:
            print(f"ERROR updating {f.name}: {'; '.join(errs)}", flush=True)
            sys.exit(1)

    else:
        # Non-parallel mode: process all files sequentially
        print(f"Processing {total} FITS files in {BASE_DIR}", flush=True)

        n_ok = n_err = 0
        error_files = []

        # Shuffle so we can parallelize
        import random
        random.shuffle(fits_list)

        for ii, f in enumerate(fits_list):
            print(f"Processing {f} ({ii} of {len(fits_list)})", flush=True)
            ok, errs = apply_fixes(f)
            if ok:
                n_ok += 1
            else:
                n_err += 1
                error_files.append((f.name, errs))

        print("\n--- Summary ---", flush=True)
        print(f"Total files: {total}", flush=True)
        print(f"Successfully updated: {n_ok}", flush=True)
        print(f"Errors: {n_err}", flush=True)

        if error_files:
            print("\nFiles with errors:", flush=True)
            for name, errs in error_files:
                print(f"  {name}: {'; '.join(errs)}", flush=True)


if __name__ == "__main__":
    main()