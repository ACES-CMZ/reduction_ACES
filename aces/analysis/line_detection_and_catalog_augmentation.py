"""
Identify spectral line detections toward ACES compact continuum sources
and augment the catalog with fitted line parameters.

For each source with extracted spectra, this script:
1. Loads the extracted 1D spectrum for each SPW
2. Checks for candidate line detections near rest frequencies from the linelist
3. Fits Gaussian profiles to promising candidates using astropy modeling
4. Optionally extracts a background annulus spectrum from the parent cube
   to confirm the detection is source-associated (not ambient emission)
5. Records fit parameters (peak, width, centroid) for common lines
   as new catalog columns, and rare detections in a JSON file.
"""

import glob
import hashlib
import json
import os
import pickle
import re
import warnings

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.io import fits
from astropy.modeling import fitting, models
from astropy.stats import mad_std, sigma_clip
from astropy.table import Table
from astropy.wcs import WCS
from astropy import log

from aces import conf

basepath = conf.basepath


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
CATALOG_NAME_PREFIX = 'ACEScatalog_v0_20260130'
SPECTRUM_DIR = os.path.join(basepath, 'spectra')
CATALOG_PATH = os.path.join(
    basepath,
    'upload/compact_cont_source_catalog/'
    'official_cont_catalog_files/aces_compact_catalog_v0.fits'
)
LINELIST_PATH = os.path.join(
    basepath,
    'reduction_ACES/aces/data/tables/linelist.csv'
)
OUTPUT_CATALOG_PATH = os.path.join(
    basepath,
    'catalogs/aces_compact_catalog_v0_withLines.fits'
)
OUTPUT_CATALOG_ECSV_PATH = os.path.join(
    basepath,
    'catalogs/aces_compact_catalog_v0_withLines.ecsv'
)
SPARSE_DETECTIONS_JSON = os.path.join(
    basepath,
    'catalogs/aces_compact_catalog_v0_rare_line_detections.json'
)

# Cache directory for intermediate results
CACHE_DIR = os.path.join(
    basepath,
    'reduction_ACES/aces/analysis/.cache'
)
CACHE_DETECTIONS = os.path.join(CACHE_DIR, 'all_detections.pkl')
CACHE_LINE_COUNTS = os.path.join(CACHE_DIR, 'line_counts.pkl')
CACHE_METADATA = os.path.join(CACHE_DIR, 'cache_metadata.json')

# Velocity range for the Galactic center (km/s)
VELO_RANGE_KMS = 200.0

# S/N threshold for a candidate detection
CANDIDATE_SNR = 5.0

# S/N threshold for source-vs-background confirmation
SOURCE_VS_BG_FACTOR = 2.0

# Minimum number of sources with a detection for a line to get a catalog
# column; otherwise it goes into the JSON sparse store
MIN_DETECTIONS_FOR_COLUMN = 5

# SPW -> 12m SPW mapping (from linelist)
LINE_SPWS = {'spw25': 25, 'spw27': 27, 'spw29': 29,
             'spw31': 31, 'spw33': 33, 'spw35': 35}


# --------------------------------------------------------------------------
# Utility functions
# --------------------------------------------------------------------------

def freq_to_velocity(freq, rest_freq):
    """Convert frequency to radio velocity (km/s) relative to rest_freq."""
    return (const.c * (1.0 - freq / rest_freq)).to(u.km / u.s).value


def velocity_to_freq(vel_kms, rest_freq_ghz):
    """Convert radio velocity (km/s) to frequency (GHz)."""
    return rest_freq_ghz * (1.0 - vel_kms / const.c.to(u.km / u.s).value)


def load_spectrum(filepath):
    """Load a 1D extracted spectrum FITS file.

    Returns
    -------
    freq_ghz : array
        Frequency axis in GHz.
    flux : array
        Flux values (Jy/beam or Jy).
    header : fits.Header
        FITS header.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            hdu = fits.open(filepath)[0]
            header = hdu.header
            data = hdu.data.squeeze()  # collapse degenerate spatial axes

            ww = WCS(header).spectral
            naxis = data.shape[0]
            pixel_indices = np.arange(naxis)
            freq_hz = ww.pixel_to_world(pixel_indices).to(u.Hz).value
            freq_ghz = freq_hz / 1.0e9

        return freq_ghz, data, header
    except Exception as ex:
        # need to print the bad filename so we can find the dang thing... but always re-raise!
        print(f"Error loading spectrum from {filepath}: {ex}", flush=True)
        raise ex


def estimate_rms(flux, line_channels_mask=None):
    """Estimate the RMS noise from line-free channels.

    Parameters
    ----------
    flux : array
        1D flux array.
    line_channels_mask : array of bool, optional
        True where line emission is expected. These channels are excluded.

    Returns
    -------
    rms : float
    """
    if line_channels_mask is not None:
        line_free = flux[~line_channels_mask]
    else:
        line_free = flux

    # Iterative sigma clipping: reject >3sigma
    for _ in range(3):
        med = np.nanmedian(line_free)
        std = np.nanstd(line_free)
        if std == 0:
            return 0.0
        good = np.abs(line_free - med) < 3.0 * std
        line_free = line_free[good]

    return np.nanstd(line_free)


def find_channels_near_line(freq_ghz, rest_freq_ghz, velo_range_kms):
    """Return a boolean mask of channels within +/-velo_range_kms of rest_freq."""
    vel_kms = freq_to_velocity(freq_ghz * u.GHz,
                               rest_freq_ghz * u.GHz)
    return np.abs(vel_kms) <= velo_range_kms


def fit_gaussian_to_line(freq_ghz, flux, rest_freq_ghz,
                         velo_range_kms=VELO_RANGE_KMS):
    """Attempt a Gaussian + constant baseline fit near a spectral line.

    Parameters
    ----------
    freq_ghz : array
        Full frequency axis.
    flux : array
        Full flux array.
    rest_freq_ghz : float
        Rest frequency of the line in GHz.
    velo_range_kms : float
        Velocity window around the line in km/s.

    Returns
    -------
    result : dict or None
        Dictionary with keys 'amplitude', 'mean_freq_ghz', 'stddev_freq_ghz',
        'centroid_vel_kms', 'fwhm_kms', 'peak_flux', 'snr', 'baseline',
        or None if the fit fails or is not significant.
    """
    mask = find_channels_near_line(freq_ghz, rest_freq_ghz, velo_range_kms)
    if mask.sum() < 5:
        return None

    xdata = freq_ghz[mask]
    ydata = flux[mask]

    # Estimate baseline and rms using sigma clipping and MAD
    # Use astropy's sigma_clip to reject outliers (likely line emission)
    clipped = sigma_clip(ydata, sigma=3.0, maxiters=3, masked=True, stdfunc=mad_std)
    clipped_data = clipped.compressed()  # Get non-masked values
    
    if len(clipped_data) < 3:
        log.error("Not enough data points after clipping for baseline fit.")
        return None
    
    # Final baseline and RMS estimates
    baseline_est = np.nanmedian(clipped_data)
    rms = mad_std(clipped_data, ignore_nan=True)
    if rms == 0:
        log.error("RMS is zero after baseline estimation.")
        return None

    ydata_sub = ydata - baseline_est
    peak_val = np.nanmax(ydata_sub)
    peak_idx = np.nanargmax(ydata_sub)
    snr = peak_val / rms

    if snr < CANDIDATE_SNR:
        return None

    # Initial guess for Gaussian
    amp_init = peak_val
    mean_init = xdata[peak_idx]
    # Estimate width from half-maximum
    above_half = ydata_sub > 0.5 * peak_val
    if above_half.sum() > 1:
        channel_width = np.abs(np.diff(xdata).mean())
        stddev_init = above_half.sum() * channel_width / 2.35
    else:
        stddev_init = np.abs(np.diff(xdata).mean()) * 3

    # Check for non-finite values in input data
    finite_mask = np.isfinite(xdata) & np.isfinite(ydata)
    if not np.any(finite_mask):
        log.error("No finite data points available for fitting.")
        return None
    
    # Use only finite data
    xdata_clean = xdata[finite_mask]
    ydata_clean = ydata[finite_mask]
    
    gauss = models.Gaussian1D(amplitude=amp_init, mean=mean_init,
                              stddev=stddev_init)
    fitter = fitting.LevMarLSQFitter()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fitted = fitter(gauss, xdata_clean, ydata_clean)
    fit_amp = fitted.amplitude.value
    fit_mean = fitted.mean.value
    fit_stddev = abs(fitted.stddev.value)
    # Sanity checks on fit
    # Mean must be within the velocity window
    fit_vel = freq_to_velocity(fit_mean * u.GHz,
                               rest_freq_ghz * u.GHz)
    if abs(fit_vel) > velo_range_kms:
        return None
    # Width must be positive and not absurdly large
    fwhm_ghz = 2.3548 * fit_stddev
    fwhm_vel = abs(freq_to_velocity((rest_freq_ghz + fwhm_ghz / 2) * u.GHz,
                                    rest_freq_ghz * u.GHz) -
                   freq_to_velocity((rest_freq_ghz - fwhm_ghz / 2) * u.GHz,
                                    rest_freq_ghz * u.GHz))
    if fwhm_vel > 300 or fwhm_vel < 0.5:
        return None

    # Re-estimate S/N of the fit
    residuals = ydata - fitted(xdata)
    rms_resid = np.nanstd(residuals)
    if rms_resid == 0:
        log.error("RMS of residuals is zero after fitting.")
        fit_snr = snr
    else:
        fit_snr = fit_amp / rms_resid

    if fit_snr < CANDIDATE_SNR:
        return None

    return {
        'amplitude': fit_amp,
        'mean_freq_ghz': fit_mean,
        'stddev_freq_ghz': fit_stddev,
        'centroid_vel_kms': fit_vel,
        'fwhm_kms': fwhm_vel,
        'peak_flux': fit_amp + baseline_est,
        'snr': fit_snr,
        'baseline': baseline_est,
    }


def load_linelist(filepath=LINELIST_PATH):
    """Load the ACES linelist CSV and return a table with parsed line info.

    Returns
    -------
    linelist : astropy Table
        Columns include: 'Line', 'Rest (GHz)', '12m SPW', 'col9' (importance)
    """
    tbl = Table.read(filepath, format='ascii.csv')

    # Clean up the rest frequency column
    rest_freqs = []
    for row in tbl:
        rf = row['Rest (GHz)']
        if isinstance(rf, str):
            rf = rf.strip()
        rest_freqs.append(float(rf))
    tbl['rest_freq_ghz'] = rest_freqs

    # Create a clean line name for column naming
    clean_names = []
    for row in tbl:
        name = row['Line'].strip()
        # Replace problematic characters for column names
        cname = re.sub(r'[^a-zA-Z0-9_]', '_', name)
        cname = re.sub(r'_+', '_', cname).strip('_')
        clean_names.append(cname)
    tbl['clean_name'] = clean_names

    return tbl


def get_source_spectra_files(source_id, spectrum_dir=SPECTRUM_DIR,
                             prefix=CATALOG_NAME_PREFIX):
    """Find all extracted spectrum FITS files for a given source ID.

    Returns
    -------
    spw_files : dict
        {spw_label: filepath} e.g. {'spw25': '/path/to/file.fits'}
    """
    pattern = os.path.join(
        spectrum_dir,
        f"{prefix}_source{source_id}_ellipseaverage_*.fits"
    )
    files = glob.glob(pattern)
    spw_files = {}
    for fn in files:
        for spw_label in LINE_SPWS:
            if f'.{spw_label}.' in fn:
                spw_files[spw_label] = fn
                break
    return spw_files


def extract_background_spectrum_from_cube(cubefn, center_glon, center_glat,
                                          inner_radius_arcsec,
                                          outer_radius_arcsec):
    """Extract a mean spectrum from an annulus around a source position.

    This extracts directly from the cube to get the background emission level.

    Parameters
    ----------
    cubefn : str
        Path to the spectral cube FITS file.
    center_glon, center_glat : float
        Galactic coordinates of the source center (degrees).
    inner_radius_arcsec : float
        Inner radius of the annulus (arcsec).
    outer_radius_arcsec : float
        Outer radius of the annulus (arcsec).

    Returns
    -------
    freq_ghz : array
        Frequency axis in GHz.
    bg_spectrum : array
        Mean background spectrum.
    """
    import regions
    from astropy.coordinates import SkyCoord
    from spectral_cube import SpectralCube

    center = SkyCoord(center_glon, center_glat,
                      frame='galactic', unit=(u.deg, u.deg))

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cube = SpectralCube.read(cubefn)
        cube.allow_huge_operations = True

    # Use inner and outer circular apertures to make an annulus
    inner_reg = regions.CircleSkyRegion(
        center, radius=inner_radius_arcsec * u.arcsec)
    outer_reg = regions.CircleSkyRegion(
        center, radius=outer_radius_arcsec * u.arcsec)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        outer_scube = cube.subcube_from_regions([outer_reg])
        outer_mean = outer_scube.mean(axis=(1, 2))

        inner_scube = cube.subcube_from_regions([inner_reg])
        inner_mean = inner_scube.mean(axis=(1, 2))

    # Annulus = (outer * N_outer - inner * N_inner) / (N_outer - N_inner)
    # For simplicity, approximate as outer average when inner is small
    freq_hz = cube.spectral_axis.to(u.Hz).value
    freq_ghz = freq_hz / 1e9

    # Get pixel counts for proper annulus weighting
    outer_mask = outer_scube.mask.include()
    inner_mask = inner_scube.mask.include()
    n_outer = np.sum(outer_mask, axis=(1, 2)).astype(float)
    n_inner = np.sum(inner_mask, axis=(1, 2)).astype(float)
    n_annulus = n_outer - n_inner
    n_annulus[n_annulus <= 0] = 1.0

    bg_spectrum = ((outer_mean.value * n_outer - inner_mean.value * n_inner)
                   / n_annulus)

    return freq_ghz, bg_spectrum


def find_parent_cube(source_id, spw_label, spectrum_dir=SPECTRUM_DIR,
                     prefix=CATALOG_NAME_PREFIX):
    """Find the parent cube path from the spectrum filename.

    The spectrum filename encodes the cube UID, so we reconstruct the cube path.
    """
    spw_files = get_source_spectra_files(source_id, spectrum_dir, prefix)
    if spw_label not in spw_files:
        return None

    spec_fn = spw_files[spw_label]
    basename = os.path.basename(spec_fn)
    # Strip prefix and source info to get the cube filename
    # Format: PREFIX_sourceNNN_ellipseaverage_CUBENAME
    parts = basename.split('_ellipseaverage_')
    if len(parts) < 2:
        return None
    cube_basename = parts[1]

    product_dir = os.path.join(
        basepath,
        'rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/'
        'group.uid___A001_X1590_X30a9/'
    )

    # Extract the MOUS uid from the cube filename
    # e.g. uid___A001_X15a0_Xf4.s38_0.Sgr_A_star_sci.spw29.cube...
    uid_match = re.match(r'(uid___A001_X\w+_X\w+)', cube_basename)
    if uid_match is None:
        return None
    uid = uid_match.group(1)
    mous_part = uid.replace('uid___A001_', '')

    # Search for the cube
    pattern = os.path.join(
        product_dir,
        f'member.uid___A001_{mous_part}',
        'calibrated/working',
        cube_basename.replace('.pbcor.fits', '.pbcor.fits')
    )
    candidates = glob.glob(pattern)
    if candidates:
        return candidates[0]

    # Try broader search
    pattern2 = os.path.join(
        product_dir,
        f'member.uid___A001_{mous_part}',
        'calibrated/working',
        f'*{spw_label}*pbcor.fits'
    )
    candidates2 = glob.glob(pattern2)
    if candidates2:
        return candidates2[0]

    return None


# --------------------------------------------------------------------------
# Cache utilities
# --------------------------------------------------------------------------

def compute_file_hash(filepath):
    """Compute MD5 hash of a file for cache invalidation.
    
    Parameters
    ----------
    filepath : str
        Path to file.
        
    Returns
    -------
    hash_str : str
        Hexadecimal hash string.
    """
    hasher = hashlib.md5()
    with open(filepath, 'rb') as fh:
        # Read in chunks to handle large files
        for chunk in iter(lambda: fh.read(8192), b''):
            hasher.update(chunk)
    return hasher.hexdigest()


def get_cache_metadata():
    """Load cache metadata or return empty dict if not found."""
    if os.path.exists(CACHE_METADATA):
        with open(CACHE_METADATA, 'r') as fh:
            return json.load(fh)
    return {}


def save_cache_metadata(metadata):
    """Save cache metadata to disk."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(CACHE_METADATA, 'w') as fh:
        json.dump(metadata, fh, indent=2)


def is_cache_valid():
    """Check if cached results are valid based on input file hashes.
    
    Returns
    -------
    valid : bool
        True if cache exists and input files haven't changed.
    """
    if not os.path.exists(CACHE_DETECTIONS):
        return False
    
    metadata = get_cache_metadata()
    if not metadata:
        return False
    
    # Check if input files have changed
    current_catalog_hash = compute_file_hash(CATALOG_PATH)
    current_linelist_hash = compute_file_hash(LINELIST_PATH)
    
    if metadata.get('catalog_hash') != current_catalog_hash:
        print("  Cache invalid: catalog file has changed", flush=True)
        return False
    
    if metadata.get('linelist_hash') != current_linelist_hash:
        print("  Cache invalid: linelist file has changed", flush=True)
        return False
    
    # Check if this script has been modified
    script_path = __file__
    current_script_hash = compute_file_hash(script_path)
    if metadata.get('script_hash') != current_script_hash:
        print("  Cache invalid: analysis script has changed", flush=True)
        return False
    
    return True


def save_detections_cache(all_detections, line_counts):
    """Save detection results to cache.
    
    Parameters
    ----------
    all_detections : dict
        Detection results.
    line_counts : dict
        Line detection counts.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    
    print(f"\nSaving detections cache to {CACHE_DETECTIONS}", flush=True)
    with open(CACHE_DETECTIONS, 'wb') as fh:
        pickle.dump(all_detections, fh)
    
    print(f"Saving line counts cache to {CACHE_LINE_COUNTS}", flush=True)
    with open(CACHE_LINE_COUNTS, 'wb') as fh:
        pickle.dump(line_counts, fh)
    
    # Update metadata
    metadata = {
        'catalog_hash': compute_file_hash(CATALOG_PATH),
        'linelist_hash': compute_file_hash(LINELIST_PATH),
        'script_hash': compute_file_hash(__file__),
        'n_sources_with_detections': len(all_detections),
        'total_detections': sum(len(d) for d in all_detections.values()),
    }
    save_cache_metadata(metadata)
    print("Cache metadata saved", flush=True)


def load_detections_cache():
    """Load detection results from cache.
    
    Returns
    -------
    all_detections : dict
        Detection results.
    line_counts : dict
        Line detection counts.
    """
    print(f"\nLoading detections from cache: {CACHE_DETECTIONS}", flush=True)
    with open(CACHE_DETECTIONS, 'rb') as fh:
        all_detections = pickle.load(fh)
    
    print(f"Loading line counts from cache: {CACHE_LINE_COUNTS}", flush=True)
    with open(CACHE_LINE_COUNTS, 'rb') as fh:
        line_counts = pickle.load(fh)
    
    return all_detections, line_counts


# --------------------------------------------------------------------------
# Main analysis
# --------------------------------------------------------------------------

def scan_all_sources(catalog, linelist, spectrum_dir=SPECTRUM_DIR,
                     prefix=CATALOG_NAME_PREFIX):
    """Scan all catalog sources for spectral line detections.

    Parameters
    ----------
    catalog : astropy Table
        The ACES compact source catalog.
    linelist : astropy Table
        The linelist with rest frequencies.

    Returns
    -------
    all_detections : dict
        {source_index: {line_clean_name: fit_result_dict}}
    """
    # Group linelist by SPW for efficient processing
    lines_by_spw = {}
    for ll_row in linelist:
        spw_num = ll_row['12m SPW']
        spw_label = f'spw{spw_num}'
        if spw_label not in lines_by_spw:
            lines_by_spw[spw_label] = []
        lines_by_spw[spw_label].append(ll_row)
    
    all_detections = {}
    n_sources = len(catalog)

    for ii, row in enumerate(catalog):
        src_id = row['index']
        spw_files = get_source_spectra_files(src_id, spectrum_dir, prefix)

        if not spw_files:
            continue

        source_detections = {}

        # Loop over each SPW that has a spectrum file
        for spw_label, spw_filepath in spw_files.items():
            # Load the spectrum once for this SPW
            freq_ghz, flux, header = load_spectrum(spw_filepath)
            
            # Check if we have any lines in this SPW
            if spw_label not in lines_by_spw:
                continue
            
            # Process all lines in this SPW
            for ll_row in lines_by_spw[spw_label]:
                rest_freq = ll_row['rest_freq_ghz']
                
                # Skip if the rest frequency is not within the spectral range
                if rest_freq < freq_ghz.min() or rest_freq > freq_ghz.max():
                    continue

                result = fit_gaussian_to_line(freq_ghz, flux, rest_freq)
                if result is not None:
                    clean_name = ll_row['clean_name']
                    source_detections[clean_name] = result
                    source_detections[clean_name]['line_name'] = ll_row['Line']
                    source_detections[clean_name]['rest_freq_ghz'] = rest_freq
                    source_detections[clean_name]['spw'] = spw_label

        if source_detections:
            all_detections[int(src_id)] = source_detections

        if (ii + 1) % 100 == 0:
            n_det = len(all_detections)
            print(f"Processed {ii + 1}/{n_sources} sources, "
                  f"{n_det} with detections so far", flush=True)

    return all_detections


def verify_source_detections(all_detections, catalog, linelist,
                             spectrum_dir=SPECTRUM_DIR,
                             prefix=CATALOG_NAME_PREFIX):
    """For sources with promising detections, check whether the line emission
    is genuinely source-associated by comparing to a background annulus.

    We re-extract a background annulus spectrum from the parent cube and
    compare the line peak in the source aperture to the background.

    Parameters
    ----------
    all_detections : dict
        Output of scan_all_sources.
    catalog : astropy Table
        The ACES catalog.
    linelist : astropy Table
        The linelist.

    Returns
    -------
    verified_detections : dict
        Same structure as all_detections but only including lines whose
        source peak exceeds background peak by SOURCE_VS_BG_FACTOR.
    """
    verified = {}

    catalog_by_idx = {row['index']: row for row in catalog}

    for src_id, detections in all_detections.items():
        print(src_id, end=', ')
        row = catalog_by_idx[src_id]
        source_major = row['fitted_major']  # arcsec
        source_minor = row['fitted_minor']  # arcsec
        # Use larger of the two axes as inner radius, 3x as outer
        inner_r = max(source_major, source_minor) * 1.5
        outer_r = inner_r * 3.0

        verified_lines = {}
        for line_name, det in detections.items():
            spw_label = det['spw']
            cubefn = find_parent_cube(src_id, spw_label,
                                      spectrum_dir, prefix)
            if cubefn is None:
                # Cannot verify without the cube; keep it but flag
                det['verified'] = False
                verified_lines[line_name] = det
                continue

            bg_freq, bg_spec = extract_background_spectrum_from_cube(
                cubefn,
                row['GLON_peak'], row['GLAT_peak'],
                inner_r, outer_r
            )

            # Find the peak in background at the same frequency range
            rest_freq = det['rest_freq_ghz']
            bg_mask = find_channels_near_line(
                bg_freq, rest_freq, VELO_RANGE_KMS)
            if bg_mask.sum() < 3:
                det['verified'] = False
                verified_lines[line_name] = det
                continue

            bg_baseline = np.nanmedian(bg_spec[~bg_mask]) if (~bg_mask).sum() > 0 else 0
            bg_peak = np.nanmax(bg_spec[bg_mask]) - bg_baseline

            source_excess = det['amplitude']
            if bg_peak <= 0:
                ratio = np.inf
            else:
                ratio = source_excess / bg_peak

            det['bg_peak'] = float(bg_peak)
            det['source_to_bg_ratio'] = float(ratio)
            det['verified'] = ratio >= SOURCE_VS_BG_FACTOR

            if det['verified']:
                verified_lines[line_name] = det

        if verified_lines:
            verified[src_id] = verified_lines

    print()
    return verified


def tally_detections(all_detections):
    """Count how many sources have a detection for each line.

    Returns
    -------
    line_counts : dict
        {clean_line_name: count}
    """
    counts = {}
    for src_id, dets in all_detections.items():
        for line_name in dets:
            counts[line_name] = counts.get(line_name, 0) + 1
    return counts


def augment_catalog(catalog, all_detections, line_counts):
    """Add line detection columns to the catalog.

    Lines detected in >= MIN_DETECTIONS_FOR_COLUMN sources get three columns
    each: {line}_peak, {line}_fwhm_kms, {line}_centroid_kms.

    Rare detections are returned as a dict to be written to JSON.

    Parameters
    ----------
    catalog : astropy Table
        Original catalog.
    all_detections : dict
        Detection results.
    line_counts : dict
        {line_name: count}.

    Returns
    -------
    catalog : astropy Table
        Augmented catalog (modified in-place).
    sparse_detections : dict
        {source_id: {line_name: {peak, fwhm_kms, centroid_kms, ...}}}
    """
    common_lines = sorted(
        [ln for ln, cnt in line_counts.items()
         if cnt >= MIN_DETECTIONS_FOR_COLUMN]
    )
    rare_lines = sorted(
        [ln for ln, cnt in line_counts.items()
         if cnt < MIN_DETECTIONS_FOR_COLUMN]
    )

    print(f"\nCommon lines (>= {MIN_DETECTIONS_FOR_COLUMN} detections, "
          f"get catalog columns): {len(common_lines)}")
    for ln in common_lines:
        print(f"  {ln}: {line_counts[ln]} detections")

    print(f"\nRare lines (< {MIN_DETECTIONS_FOR_COLUMN} detections, "
          f"stored in JSON): {len(rare_lines)}")
    for ln in rare_lines:
        print(f"  {ln}: {line_counts[ln]} detections")

    # Initialize columns for common lines
    n = len(catalog)
    for ln in common_lines:
        catalog[f'{ln}_peak'] = np.full(n, np.nan)
        catalog[f'{ln}_fwhm_kms'] = np.full(n, np.nan)
        catalog[f'{ln}_centroid_kms'] = np.full(n, np.nan)
        catalog[f'{ln}_snr'] = np.full(n, np.nan)

    # Build an index map: catalog index -> row position
    idx_to_row = {}
    for ii, row in enumerate(catalog):
        idx_to_row[row['index']] = ii

    sparse_detections = {}

    for src_id, dets in all_detections.items():
        if src_id not in idx_to_row:
            continue
        row_idx = idx_to_row[src_id]

        for line_name, det in dets.items():
            if line_name in common_lines:
                catalog[f'{line_name}_peak'][row_idx] = det['amplitude']
                catalog[f'{line_name}_fwhm_kms'][row_idx] = det['fwhm_kms']
                catalog[f'{line_name}_centroid_kms'][row_idx] = det['centroid_vel_kms']
                catalog[f'{line_name}_snr'][row_idx] = det['snr']
            elif line_name in rare_lines:
                src_key = str(src_id)
                if src_key not in sparse_detections:
                    sparse_detections[src_key] = {}
                sparse_detections[src_key][line_name] = {
                    'peak': float(det['amplitude']),
                    'fwhm_kms': float(det['fwhm_kms']),
                    'centroid_kms': float(det['centroid_vel_kms']),
                    'snr': float(det['snr']),
                    'line_name': det.get('line_name', line_name),
                    'rest_freq_ghz': float(det['rest_freq_ghz']),
                }

    return catalog, sparse_detections


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

if __name__ == "__main__":
    print("Loading catalog...", flush=True)
    catalog = Table.read(CATALOG_PATH)
    print(f"  {len(catalog)} sources", flush=True)

    print("Loading linelist...", flush=True)
    linelist = load_linelist()
    print(f"  {len(linelist)} lines across {len(set(linelist['12m SPW']))} SPWs",
          flush=True)

    # Filter linelist to only lines with importance markers (**, ***) for
    # the initial scan to save time; the full list has 72 lines
    # We scan all lines but prioritize the important ones
    print("\nLines by importance:")
    for imp in ('***', '**', '*', ''):
        mask = linelist['col9'] == imp
        if mask.sum() > 0:
            label = imp if imp else '(unmarked)'
            print(f"  {label}: {mask.sum()} lines")

    # Check if we can use cached results
    print("\nChecking for cached detection results...", flush=True)
    if is_cache_valid():
        print("  Cache is valid, loading cached results...", flush=True)
        all_detections, line_counts = load_detections_cache()
        n_sources_with_det = len(all_detections)
        total_dets = sum(len(d) for d in all_detections.values())
        print(f"  Loaded: {n_sources_with_det} sources with "
              f"{total_dets} total line detections", flush=True)
    else:
        print("  No valid cache found, scanning all sources...", flush=True)
        print("\nScanning all sources for line detections...", flush=True)
        all_detections = scan_all_sources(catalog, linelist)
        n_sources_with_det = len(all_detections)
        total_dets = sum(len(d) for d in all_detections.values())
        print(f"\nInitial scan complete: {n_sources_with_det} sources with "
              f"{total_dets} total line detections", flush=True)

        # Tally and save to cache
        line_counts = tally_detections(all_detections)
        save_detections_cache(all_detections, line_counts)
    print("\nDetection counts per line:")
    for ln, cnt in sorted(line_counts.items(), key=lambda x: -x[1]):
        print(f"  {ln}: {cnt}")

    # Augment catalog
    print("\nAugmenting catalog...", flush=True)
    catalog, sparse_dets = augment_catalog(catalog, all_detections,
                                           line_counts)

    # Write outputs
    print(f"\nWriting augmented catalog to {OUTPUT_CATALOG_PATH}", flush=True)
    catalog.write(OUTPUT_CATALOG_PATH, overwrite=True)
    print(f"Writing augmented catalog to {OUTPUT_CATALOG_ECSV_PATH}",
          flush=True)
    catalog.write(OUTPUT_CATALOG_ECSV_PATH, format='ascii.ecsv',
                  overwrite=True)

    if sparse_dets:
        print(f"Writing {len(sparse_dets)} rare-line source entries to "
              f"{SPARSE_DETECTIONS_JSON}", flush=True)
        with open(SPARSE_DETECTIONS_JSON, 'w') as fh:
            json.dump(sparse_dets, fh, indent=2)

    print("\nDone.", flush=True)
