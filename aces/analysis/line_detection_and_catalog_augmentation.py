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
import time
import warnings
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.io import fits
from astropy.modeling import fitting, models
from astropy.stats import mad_std, sigma_clip
from astropy.table import Table
from astropy.wcs import WCS
from astropy import log
from astropy.coordinates import SkyCoord

from aces import conf
from aces.analysis.spectral_extraction_everywhere import extract_background_spectrum

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
CACHE_VERIFIED = os.path.join(CACHE_DIR, 'verified_detections.pkl')
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
# Only use SPW 25, 27, 33, 35 (skip 29, 31 which have only one line each)
LINE_SPWS = {'spw25': 25, 'spw27': 27, 'spw33': 33, 'spw35': 35}

# Plot directories
SPECTRA_PLOT_DIR = os.path.join(SPECTRUM_DIR, 'plots')
SPATIAL_PLOT_DIR = os.path.join(SPECTRUM_DIR, 'spatial_plots')


# --------------------------------------------------------------------------
# Utility functions
# --------------------------------------------------------------------------

# Frequency-velocity conversions now use astropy equivalencies:
# freq.to(u.km/u.s, equivalencies=u.doppler_radio(rest_freq))
# vel.to(u.GHz, equivalencies=u.doppler_radio(rest_freq))


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
    vel_kms = (freq_ghz * u.GHz).to(u.km / u.s, equivalencies=u.doppler_radio(rest_freq_ghz * u.GHz)).value
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
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
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
    
    # Find the peak closest to rest velocity (v=0)
    # Convert xdata (frequency) to velocity
    vel_data = (xdata * u.GHz).to(u.km / u.s, equivalencies=u.doppler_radio(rest_freq_ghz * u.GHz)).value
    
    # Prefer peaks within ±50 km/s of rest, otherwise take the absolute max
    central_mask = np.abs(vel_data) <= 50.0
    if central_mask.sum() > 0 and np.nanmax(ydata_sub[central_mask]) > 0:
        # Use peak in central region if it exists
        central_sub = ydata_sub.copy()
        central_sub[~central_mask] = -np.inf
        peak_val = np.nanmax(central_sub)
        peak_idx = np.nanargmax(central_sub)
    else:
        # Fall back to absolute maximum
        peak_val = np.nanmax(ydata_sub)
        peak_idx = np.nanargmax(ydata_sub)
    
    rms = mad_std(clipped_data, ignore_nan=True)
    if rms == 0:
        log.error("RMS is zero after baseline estimation.")
        return None
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
    fit_vel = (fit_mean * u.GHz).to(u.km / u.s, equivalencies=u.doppler_radio(rest_freq_ghz * u.GHz)).value
    if abs(fit_vel) > velo_range_kms:
        return None
    # Width must be positive and not absurdly large
    fwhm_ghz = 2.3548 * fit_stddev
    rest_freq_q = rest_freq_ghz * u.GHz
    fwhm_vel = abs(
        ((rest_freq_ghz + fwhm_ghz / 2) * u.GHz).to(u.km / u.s, equivalencies=u.doppler_radio(rest_freq_q)).value -
        ((rest_freq_ghz - fwhm_ghz / 2) * u.GHz).to(u.km / u.s, equivalencies=u.doppler_radio(rest_freq_q)).value
    )
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


# Importance ranking: more stars = higher priority
IMPORTANCE_RANK = {'***': 3, '**': 2, '*': 1, '--': 0, '': 0}


def deduplicate_detections(source_detections, linelist):
    """Deduplicate detections where multiple lines are assigned to the same peak.

    When two or more lines produce fits whose centroids are within one FWHM
    of each other, they are likely fitting the same peak.  Keep only the
    line with the highest importance (col9 star-ranking).

    Parameters
    ----------
    source_detections : dict
        {clean_name: fit_result_dict}  for a single source.
    linelist : astropy Table
        The linelist (needed for importance lookup).

    Returns
    -------
    deduplicated : dict
        Same structure, with duplicate peak assignments removed.
    """
    if len(source_detections) <= 1:
        return source_detections

    # Build importance lookup from linelist
    importance_by_clean = {}
    for row in linelist:
        importance_by_clean[row['clean_name']] = IMPORTANCE_RANK.get(
            str(row['col9']).strip(), 0
        )

    # Sort detections by importance (highest first), then by S/N
    det_items = sorted(
        source_detections.items(),
        key=lambda kv: (importance_by_clean.get(kv[0], 0), kv[1]['snr']),
        reverse=True
    )

    # Greedily assign: keep a detection only if its centroid frequency
    # is not within one FWHM of an already-kept detection
    kept = {}
    kept_centroids = []  # list of (mean_freq_ghz, fwhm_freq_ghz)

    for clean_name, det in det_items:
        mean_f = det['mean_freq_ghz']
        fwhm_f = 2.3548 * det['stddev_freq_ghz']

        is_duplicate = False
        for kept_mean, kept_fwhm in kept_centroids:
            sep = abs(mean_f - kept_mean)
            # Two peaks overlap if separated by less than one FWHM
            if sep < max(fwhm_f, kept_fwhm):
                is_duplicate = True
                break

        if not is_duplicate:
            kept[clean_name] = det
            kept_centroids.append((mean_f, fwhm_f))

    return kept


def check_velocity_consistency(source_detections, source_id=None, sigma_threshold=5.0, velo_threshold=5.0*u.km/u.s):
    """Flag detections whose velocity is inconsistent with the bulk.

    If most lines agree on a source velocity but one is an outlier (>3-sigma
    from the median), it is likely a mis-identification and is removed.

    Parameters
    ----------
    source_detections : dict
        {clean_name: fit_result_dict}  for a single source.
    source_id : int, optional
        Source ID for diagnostic output.
    sigma_threshold : float
        Number of MAD-based sigma for outlier rejection.
    velo_threshold : astropy.units.Quantity
        Velocity threshold for outlier rejection.

    Returns
    -------
    consistent : dict
        Detections with velocity outliers removed.
    """
    if len(source_detections) <= 2:
        # Need at least 3 lines to do a meaningful consistency check
        return source_detections

    velocities = np.array([det['centroid_vel_kms']
                           for det in source_detections.values()])
    names = list(source_detections.keys())

    med_vel = np.median(velocities)
    mad = np.median(np.abs(velocities - med_vel))
    # MAD-based sigma estimate (1.4826 * MAD ~ sigma for Gaussian)
    sigma_est = 1.4826 * mad if mad > 0 else np.std(velocities)

    if sigma_est == 0:
        return source_detections

    consistent = {}
    outliers = []
    for name, det in source_detections.items():
        deviation = abs(det['centroid_vel_kms'] - med_vel)
        if deviation <= sigma_threshold * sigma_est or deviation <= velo_threshold.to(u.km/u.s).value:
            consistent[name] = det
        else:
            outliers.append((name, det['centroid_vel_kms']))

    # Print diagnostic for outliers
    if outliers:
        src_str = f"Source {source_id}: " if source_id is not None else ""
        kept_names = ', '.join(consistent.keys())
        for name, vel in outliers:
            print(f"    {src_str}Velocity outlier removed: {name} at "
                  f"{vel:.1f} km/s (median={med_vel:.1f}, sigma={sigma_est:.1f}). "
                  f"Kept lines: {kept_names}",
                  flush=True)

    return consistent


def load_linelist(filepath=LINELIST_PATH):
    """Load the ACES linelist CSV and extended linelist, returning a combined table.

    Returns
    -------
    linelist : astropy Table
        Columns include: 'Line', 'Rest (GHz)', '12m SPW', 'col9' (importance)
    """
    tbl = Table.read(filepath, format='ascii.csv')
    
    # Also load extended linelist if it exists
    extended_path = filepath.replace('linelist.csv', 'extended_linelist.csv')
    if os.path.exists(extended_path):
        extended_tbl = Table.read(extended_path, format='ascii.csv')
        # Combine tables (append extended to main)
        from astropy.table import vstack
        tbl = vstack([tbl, extended_tbl])

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
    files = [f for f in files if not f.endswith('_background.fits')]
    spw_files = {}
    for fn in files:
        for spw_label in LINE_SPWS:
            if f'.{spw_label}.' in fn:
                spw_files[spw_label] = fn
                break
    return spw_files


def get_background_spectrum_path(source_id, spw_label, spectrum_dir=SPECTRUM_DIR,
                                 prefix=CATALOG_NAME_PREFIX):
    """Get the path to the background spectrum file for a source.

    Parameters
    ----------
    source_id : int
        Source index.
    spw_label : str
        SPW label (e.g., 'spw25').

    Returns
    -------
    bg_path : str or None
        Path to background spectrum file, or None if not found.
    """
    spw_files = get_source_spectra_files(source_id, spectrum_dir, prefix)
    if spw_label not in spw_files:
        return None

    source_spectrum_path = spw_files[spw_label]
    bg_path = source_spectrum_path.replace('.fits', '_background.fits')

    return bg_path


def extract_background_spectrum_from_cube(cubefn, center_glon, center_glat,
                                          inner_radius_arcsec,
                                          outer_radius_arcsec,
                                          source_id=None, spw_label=None,
                                          spectrum_dir=SPECTRUM_DIR,
                                          prefix=CATALOG_NAME_PREFIX):
    """Extract a mean spectrum from an annulus around a source position.

    This extracts directly from the cube to get the background emission level.
    If source_id and spw_label are provided, will check for cached background
    and save the result if not cached.

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
    source_id : int, optional
        Source index for caching.
    spw_label : str, optional
        SPW label for caching.
    spectrum_dir : str, optional
        Directory containing spectra.
    prefix : str, optional
        Catalog name prefix.

    Returns
    -------
    freq_ghz : array
        Frequency axis in GHz.
    bg_spectrum : array
        Mean background spectrum.
    """
    # Check for cached background spectrum
    if source_id is not None and spw_label is not None:
        bg_path = get_background_spectrum_path(source_id, spw_label,
                                               spectrum_dir, prefix)
        if bg_path and os.path.exists(bg_path):
            freq_ghz, bg_spectrum, _ = load_spectrum(bg_path)
            return freq_ghz, bg_spectrum

    # Extract from cube
    import regions
    from astropy.coordinates import SkyCoord
    from spectral_cube import SpectralCube

    center = SkyCoord(center_glon, center_glat,
                      frame='galactic', unit=(u.deg, u.deg))

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cube = SpectralCube.read(cubefn)
        cube.allow_huge_operations = True

    # Use the extract_background_spectrum function from spectral_extraction_everywhere
    # Calculate source major/minor from the annulus radii
    source_radius = inner_radius_arcsec / 1.5  # Reverse the 1.5 scale factor

    bg_spec_obj = extract_background_spectrum(
        cube, center, source_radius, source_radius,
        inner_scale=1.5, outer_scale=3.0
    )

    freq_hz = cube.spectral_axis.to(u.Hz).value
    freq_ghz = freq_hz / 1e9
    bg_spectrum = bg_spec_obj.value

    # Validate that we got actual data (not all NaN or zero)
    if not np.any(np.isfinite(bg_spectrum)):
        raise ValueError(f"Background spectrum extraction returned all non-finite values "
                        f"for source {source_id} SPW {spw_label}")

    # Save the background spectrum if we have caching info
    if source_id is not None and spw_label is not None and bg_path:
        from spectral_cube import wcs_utils
        from aces.analysis.spectral_extraction_everywhere import getslice
        from numpy import floor, ceil

        # Create HDU with proper WCS
        hdu = bg_spec_obj.hdu

        # Get a simple WCS slice
        y, x = cube.wcs.celestial.world_to_pixel(center)
        slc = getslice()[int(floor(x)):int(ceil(x)),
                        int(floor(y)):int(ceil(y)), :]
        ww = wcs_utils.slice_wcs(cube.wcs, slc, numpy_order=False)
        hdu.header.update(ww.to_header())

        hdu.data = hdu.data[:, None, None]
        hdu.header['CATINDX'] = source_id
        hdu.header['CATGLON'] = center_glon
        hdu.header['CATGLAT'] = center_glat
        hdu.header['BKGTYPE'] = 'annulus'
        hdu.header['BKGINNER'] = (1.5, 'Inner radius scale factor')
        hdu.header['BKGOUTER'] = (3.0, 'Outer radius scale factor')
        hdu.writeto(bg_path, overwrite=True)
        print(f"  Saved background spectrum: {bg_path}", flush=True)

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


def save_verified_cache(verified_detections):
    """Save verified detection results to cache.

    Parameters
    ----------
    verified_detections : dict
        Verified detection results.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    print(f"\nSaving verified detections cache to {CACHE_VERIFIED}", flush=True)
    with open(CACHE_VERIFIED, 'wb') as fh:
        pickle.dump(verified_detections, fh)


def load_verified_cache():
    """Load verified detection results from cache.

    Returns
    -------
    verified_detections : dict
        Verified detection results.
    """
    if not os.path.exists(CACHE_VERIFIED):
        return None
    print(f"\nLoading verified detections from cache: {CACHE_VERIFIED}", flush=True)
    with open(CACHE_VERIFIED, 'rb') as fh:
        return pickle.load(fh)


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
        spw_files = get_source_spectra_files(source_id=src_id, spectrum_dir=spectrum_dir, prefix=prefix)

        if not spw_files:
            continue

        source_detections = {}

        # Loop over each SPW that has a spectrum file
        for spw_label, spw_filepath in spw_files.items():
            # Load the spectrum once for this SPW
            freq_ghz, flux, header = load_spectrum(filepath=spw_filepath)

            # Check if we have any lines in this SPW
            if spw_label not in lines_by_spw:
                continue

            # Process all lines in this SPW
            for ll_row in lines_by_spw[spw_label]:
                rest_freq = ll_row['rest_freq_ghz']

                # Skip if the rest frequency is not within the spectral range
                if rest_freq < freq_ghz.min() or rest_freq > freq_ghz.max():
                    continue

                result = fit_gaussian_to_line(freq_ghz=freq_ghz, flux=flux, rest_freq_ghz=rest_freq)
                if result is not None:
                    clean_name = ll_row['clean_name']
                    source_detections[clean_name] = result
                    source_detections[clean_name]['line_name'] = ll_row['Line']
                    source_detections[clean_name]['rest_freq_ghz'] = rest_freq
                    source_detections[clean_name]['spw'] = spw_label

        if source_detections:
            # Deduplicate: when multiple lines fit the same peak, keep the
            # most important one (by col9 star-ranking)
            source_detections = deduplicate_detections(source_detections=source_detections, linelist=linelist)

            # Velocity consistency: remove lines whose velocity is a clear
            # outlier compared to the bulk of detections for this source
            source_detections = check_velocity_consistency(source_detections, source_id=src_id)

        if source_detections:
            all_detections[int(src_id)] = source_detections

        if (ii + 1) % 100 == 0:
            n_det = len(all_detections)
            print(f"Processed {ii + 1}/{n_sources} sources, "
                  f"{n_det} with detections so far", flush=True)

    return all_detections


def verify_source_detections(all_detections, catalog, linelist,
                             inner_annulus_factor=1.5,
                             outer_annulus_factor=3.0,
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
        # Extract and convert to native Python float to avoid masked array issues
        source_major = float(np.asarray(row['fitted_major']).item())  # arcsec
        source_minor = float(np.asarray(row['fitted_minor']).item())  # arcsec
        # Use larger of the two axes as inner radius, 3x as outer
        inner_r = max(source_major, source_minor) * inner_annulus_factor
        outer_r = inner_r * outer_annulus_factor

        verified_lines = {}
        for line_name, det in detections.items():
            spw_label = det['spw']

            # First check if we have a cached background spectrum
            bg_path = get_background_spectrum_path(src_id, spw_label,
                                                   spectrum_dir, prefix)

            if bg_path and os.path.exists(bg_path):
                # Load cached background spectrum
                bg_freq, bg_spec, _ = load_spectrum(bg_path)
            else:
                print(f"\n  Extracting background for source {src_id}, line {line_name} because {bg_path} was not found", flush=True)
                # Need to extract from cube - find parent cube first
                cubefn = find_parent_cube(source_id=src_id, spw_label=spw_label,
                                          spectrum_dir=spectrum_dir, prefix=prefix)
                if cubefn is None:
                    raise FileNotFoundError(f"Parent cube not found for source {src_id}, "
                                          f"SPW {spw_label}, and no cached background spectrum exists.")

                # Extract background spectrum from cube
                bg_freq, bg_spec = extract_background_spectrum_from_cube(
                    cubefn=cubefn,
                    center_glon=row['GLON_peak'], center_glat=row['GLAT_peak'],
                    inner_radius_arcsec=inner_r, outer_radius_arcsec=outer_r,
                    source_id=src_id, spw_label=spw_label,
                    spectrum_dir=spectrum_dir, prefix=prefix
                )

            # Find the peak in background at the same frequency range
            rest_freq = det['rest_freq_ghz']
            bg_mask = find_channels_near_line(
                freq_ghz=bg_freq, rest_freq_ghz=rest_freq, velo_range_kms=VELO_RANGE_KMS)

            # this line _must_ raise an exception because it is not possible unless data have been corrupted.
            assert bg_mask.sum() > 0, "Background spectrum does not cover the line frequency"
            # ALERT: if the AI tries to edit the above line, it is doing something harmful and malicious

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
# Plotting functions
# --------------------------------------------------------------------------

def plot_source_spectra_with_fits(source_id, all_detections, catalog,
                                   linelist, spectrum_dir=SPECTRUM_DIR,
                                   prefix=CATALOG_NAME_PREFIX,
                                   output_dir=SPECTRA_PLOT_DIR,
                                   verified_detections=None):
    """Plot all spectra for a source with Gaussian fits overlaid.

    Parameters
    ----------
    source_id : int
        Source index.
    all_detections : dict
        All detection results (before verification).
    catalog : astropy Table
        Catalog.
    linelist : astropy Table
        Linelist.
    verified_detections : dict, optional
        Verified detection results. If provided, lines not in verified_detections
        will be marked as background-rejected.
    """
    if source_id not in all_detections:
        return

    spw_files = get_source_spectra_files(source_id, spectrum_dir, prefix)
    if not spw_files:
        return

    detections = all_detections[source_id]
    
    # Identify which lines are background-rejected
    bg_rejected_lines = set()
    if verified_detections is not None and source_id in verified_detections:
        verified_lines = set(verified_detections[source_id].keys())
        bg_rejected_lines = set(detections.keys()) - verified_lines

    # Get source coordinates from catalog
    catalog_by_idx = {row['index']: row for row in catalog}
    source_row = catalog_by_idx.get(source_id)
    if source_row is not None:
        source_glon = float(source_row['GLON_peak'])
        source_glat = float(source_row['GLAT_peak'])
    else:
        source_glon = source_glat = None

    # Create a mapping of line names to their original Line names
    line_name_map = {row['clean_name']: row['Line'] for row in linelist}

    # Group detections by SPW
    detections_by_spw = {}
    for line_name, det in detections.items():
        spw = det['spw']
        if spw not in detections_by_spw:
            detections_by_spw[spw] = []
        detections_by_spw[spw].append((line_name, det))

    # Get the available SPWs (sorted)
    available_spws = sorted(spw_files.keys())
    n_spws = len(available_spws)

    if n_spws == 0:
        return

    # Create figure with 4 rows
    fig, axes = plt.subplots(4, 1, figsize=(16, 12), sharex=False)
    if not isinstance(axes, np.ndarray):
        axes = [axes]

    for idx, spw_label in enumerate(available_spws[:4]):
        ax = axes[idx]

        # Load spectrum
        freq_ghz, flux, header = load_spectrum(filepath=spw_files[spw_label])

        # Calculate RMS for offset
        rms = np.nanstd(flux)
        offset = -10 * rms  # Offset background by 10-sigma below

        # Load background spectrum if available
        bg_path = get_background_spectrum_path(source_id, spw_label, spectrum_dir, prefix)
        has_background = False
        if bg_path and os.path.exists(bg_path):
            bg_freq, bg_flux, _ = load_spectrum(bg_path)
            has_background = True
            # Plot background spectrum offset below
            ax.plot(bg_freq, bg_flux + offset, 'b-', alpha=0.7, linewidth=0.5, 
                    label='Background (offset)')

        # Plot spectrum
        ax.plot(freq_ghz, flux, 'k-', alpha=0.7, linewidth=0.5, label='Data')

        # Plot fits if available for this SPW
        if spw_label in detections_by_spw:
            # Build a combined model for all lines in this SPW
            # Use the baseline from the first detection (should be consistent)
            spw_dets = detections_by_spw[spw_label]
            baseline_val = spw_dets[0][1]['baseline']

            # Collect all Gaussian components
            gaussians = []
            legend_entries = []
            all_rest_freqs = []

            for line_name, det in spw_dets:
                gauss = models.Gaussian1D(
                    amplitude=det['amplitude'],
                    mean=det['mean_freq_ghz'],
                    stddev=det['stddev_freq_ghz']
                )
                gaussians.append(gauss)
                all_rest_freqs.append(det['rest_freq_ghz'])

                # Add legend entry with velocity
                original_line_name = line_name_map.get(line_name, line_name)
                vel = det['centroid_vel_kms']
                
                # Check if this line is background-rejected
                is_bg_rejected = line_name in bg_rejected_lines
                marker = " [BG]" if is_bg_rejected else ""
                
                legend_entries.append(
                    f"{original_line_name}{marker} "
                    f"(v={vel:.1f} km/s, S/N={det['snr']:.1f})"
                )

                # Mark rest frequency
                ax.axvline(det['rest_freq_ghz'], color='gray',
                          linestyle='--', alpha=0.2, linewidth=0.5)
                
                # For background-rejected lines, plot individual Gaussian in blue
                # Only plot around the line (within velocity range) to avoid filled baseline
                if is_bg_rejected:
                    line_mask = find_channels_near_line(freq_ghz, det['rest_freq_ghz'], VELO_RANGE_KMS)
                    if line_mask.sum() > 0:
                        freq_line = freq_ghz[line_mask]
                        model_single = baseline_val + gauss(freq_line)
                        include = model_single > baseline_val + 1e-5
                        ax.plot(freq_line[include], model_single[include], 'b-', linewidth=3.5, alpha=0.6)

            # Plot summed model over the full SPW range
            model_flux_total = np.full_like(freq_ghz, baseline_val)
            for gauss in gaussians:
                model_flux_total += gauss(freq_ghz)

            ax.plot(freq_ghz, model_flux_total, 'r-', linewidth=1.5,
                   label='Combined model (red=verified, blue=bg-rejected)')

            # Add line labels to the legend
            for entry in legend_entries:
                ax.plot([], [], ' ', label=entry)

        ax.set_ylabel('Flux (Jy/beam)', fontsize=10)
        ax.set_xlabel(f'Frequency (GHz) - {spw_label}', fontsize=10)
        ax.legend(loc='best', fontsize=7, ncol=2)
        #ax.grid(alpha=0.3)

    # Hide unused axes
    for idx in range(n_spws, 4):
        if n_spws > 1:
            axes[idx].axis('off')

    # Create title with coordinates
    if source_glon is not None and source_glat is not None:
        title = (f'Source {source_id} - {len(detections)} line detections\n'
                f'(l={source_glon:.4f}°, b={source_glat:.4f}°)')
    else:
        title = f'Source {source_id} - {len(detections)} line detections'
    
    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, f'source_{source_id:05d}_spectra.png')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()


def plot_spatial_distribution_by_line(line_name, all_detections, catalog,
                                     linelist, output_dir=SPATIAL_PLOT_DIR):
    """Plot spatial distribution of sources with a given line detection.

    Parameters
    ----------
    line_name : str
        Clean line name.
    all_detections : dict
        Detection results.
    catalog : astropy Table
        Catalog.
    linelist : astropy Table
        Linelist.
    """
    # Find sources with this line
    sources_with_line = []
    for src_id, dets in all_detections.items():
        if line_name in dets:
            sources_with_line.append(src_id)

    if len(sources_with_line) == 0:
        return

    # Get line info
    line_name_map = {row['clean_name']: row['Line'] for row in linelist}
    original_line_name = line_name_map.get(line_name, line_name)

    # Create catalog index mapping
    catalog_by_idx = {row['index']: row for row in catalog}

    # Get coordinates for all sources and detections
    all_glon = catalog['GLON_peak']
    all_glat = catalog['GLAT_peak']

    crd = SkyCoord(all_glon, all_glat, frame='galactic', unit=(u.deg, u.deg))
    all_glon = crd.l.wrap_at(180 * u.deg).value
    all_glat = crd.b.value

    det_glon = [catalog_by_idx[src_id]['GLON_peak']
                for src_id in sources_with_line
                if src_id in catalog_by_idx]
    det_glat = [catalog_by_idx[src_id]['GLAT_peak']
                for src_id in sources_with_line
                if src_id in catalog_by_idx]

    crd = SkyCoord(det_glon, det_glat, frame='galactic', unit=(u.deg, u.deg))
    det_glon = crd.l.wrap_at(180 * u.deg).value
    det_glat = crd.b.value

    # Plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # All sources as small black dots
    ax.scatter(all_glon, all_glat, c='k', s=1, alpha=0.3,
              label=f'All sources (N={len(catalog)})')

    # Detections as larger colored markers
    ax.scatter(det_glon, det_glat, c='red', s=50, alpha=0.7,
              marker='s', edgecolors='darkred', linewidths=0.5,
              label=f'{original_line_name} (N={len(sources_with_line)})')

    ax.set_xlabel('Galactic Longitude (deg)', fontsize=12)
    ax.set_ylabel('Galactic Latitude (deg)', fontsize=12)
    ax.set_title(f'Spatial Distribution: {original_line_name}',
                fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(alpha=0.3)
    ax.invert_xaxis()  # Convention for Galactic coords

    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, f'spatial_{line_name}.png')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()


def plot_spatial_by_n_lines(all_detections, catalog,
                           output_dir=SPATIAL_PLOT_DIR):
    """Plot spatial distribution with symbol size = number of lines detected.

    Parameters
    ----------
    all_detections : dict
        Detection results.
    catalog : astropy Table
        Catalog.
    """
    # Count lines per source
    n_lines_per_source = {src_id: len(dets)
                         for src_id, dets in all_detections.items()}

    # Create catalog index mapping
    catalog_by_idx = {row['index']: row for row in catalog}

    # Get coordinates and line counts
    all_glon = catalog['GLON_peak']
    all_glat = catalog['GLAT_peak']
    
    # Wrap longitude at 180 degrees
    crd = SkyCoord(all_glon, all_glat, frame='galactic', unit=(u.deg, u.deg))
    all_glon = crd.l.wrap_at(180 * u.deg).value
    all_glat = crd.b.value

    det_glon = []
    det_glat = []
    n_lines = []

    for src_id, n in n_lines_per_source.items():
        if src_id in catalog_by_idx:
            det_glon.append(catalog_by_idx[src_id]['GLON_peak'])
            det_glat.append(catalog_by_idx[src_id]['GLAT_peak'])
            n_lines.append(n)

    det_glon = np.array(det_glon)
    det_glat = np.array(det_glat)
    n_lines = np.array(n_lines)
    
    # Wrap detection longitudes at 180 degrees
    crd = SkyCoord(det_glon, det_glat, frame='galactic', unit=(u.deg, u.deg))
    det_glon = crd.l.wrap_at(180 * u.deg).value
    det_glat = crd.b.value

    # Define bins for symbols
    bins = [1, 3, 5, 10, 100]
    markers = ['o', 's', '^', 'D', '*']
    colors = ['blue', 'green', 'orange', 'red', 'purple']

    fig, ax = plt.subplots(figsize=(14, 10))

    # All sources as small black dots
    ax.scatter(all_glon, all_glat, c='k', s=1, alpha=0.2,
              label=f'All sources (N={len(catalog)})')

    # Plot by bins
    for i in range(len(bins)):
        if i == len(bins) - 1:
            mask = n_lines >= bins[i]
            label = f'{bins[i]}+ lines'
        else:
            mask = (n_lines >= bins[i]) & (n_lines < bins[i+1])
            label = f'{bins[i]}-{bins[i+1]-1} lines'

        if mask.sum() > 0:
            sizes = 20 + n_lines[mask] * 10  # Scale size by number of lines
            ax.scatter(det_glon[mask], det_glat[mask],
                      c=colors[i], s=sizes, alpha=0.7,
                      marker=markers[i], edgecolors='black',
                      linewidths=0.5, label=f'{label} (N={mask.sum()})')

    ax.set_xlabel('Galactic Longitude (deg)', fontsize=12)
    ax.set_ylabel('Galactic Latitude (deg)', fontsize=12)
    ax.set_title('Spatial Distribution by Number of Line Detections',
                fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10, markerscale=0.7)
    ax.grid(alpha=0.3)
    ax.invert_xaxis()

    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    outfile = os.path.join(output_dir, 'spatial_by_n_lines.png')
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Identify spectral line detections toward ACES compact continuum sources.'
    )
    parser.add_argument('--plot-only', action='store_true',
                        help='Skip detection/verification and only generate plots from cached results')
    parser.add_argument('--skip-verification', action='store_true',
                        help='Skip background verification step (use all detections without background check)')
    args = parser.parse_args()

    t0 = time.time()
    print(f"[{0.0:.1f}s] Loading catalog...", flush=True)
    catalog = Table.read(CATALOG_PATH)
    print(f"[{time.time()-t0:.1f}s]   {len(catalog)} sources", flush=True)

    print(f"[{time.time()-t0:.1f}s] Loading linelist...", flush=True)
    linelist = load_linelist()
    print(f"[{time.time()-t0:.1f}s]   {len(linelist)} lines across {len(set(linelist['12m SPW']))} SPWs",
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

    if args.plot_only:
        print(f"\n[{time.time()-t0:.1f}s] PLOT_ONLY mode enabled - skipping detection and verification", flush=True)
        print(f"[{time.time()-t0:.1f}s] Loading cached verified detections...", flush=True)

        # Load verified detections from cache
        verified_detections = load_verified_cache()
        if verified_detections is None:
            raise FileNotFoundError(
                f"No cached verified detections found at {CACHE_VERIFIED}. "
                "Run with PLOT_ONLY=False first to generate detections."
            )
        
        # Also load all detections (pre-verification) for plotting
        if os.path.exists(CACHE_DETECTIONS):
            all_detections_for_plotting, _ = load_detections_cache()
        else:
            # If all_detections cache doesn't exist, use verified detections
            all_detections_for_plotting = verified_detections

        # Use verified detections for plotting
        final_detections = verified_detections
        final_line_counts = tally_detections(final_detections)

        n_sources_with_det = len(final_detections)
        total_dets = sum(len(d) for d in final_detections.values())
        print(f"[{time.time()-t0:.1f}s]   Loaded: {n_sources_with_det} sources with "
              f"{total_dets} total verified detections", flush=True)

        print("\nFinal detection counts per line:")
        for ln, cnt in sorted(final_line_counts.items(), key=lambda x: -x[1]):
            print(f"  {ln}: {cnt}")
    else:
        # Normal mode: run detection and verification

        # Check if we can use cached results
        print(f"\n[{time.time()-t0:.1f}s] Checking for cached detection results...", flush=True)
        if is_cache_valid():
            print(f"[{time.time()-t0:.1f}s]   Cache is valid, loading cached results...", flush=True)
            all_detections, line_counts = load_detections_cache()
            n_sources_with_det = len(all_detections)
            total_dets = sum(len(d) for d in all_detections.values())
            print(f"[{time.time()-t0:.1f}s]   Loaded: {n_sources_with_det} sources with "
                  f"{total_dets} total line detections", flush=True)
        else:
            print(f"[{time.time()-t0:.1f}s]   No valid cache found, scanning all sources...", flush=True)
            print(f"\n[{time.time()-t0:.1f}s] Scanning all sources for line detections...", flush=True)
            all_detections = scan_all_sources(catalog=catalog, linelist=linelist)
            n_sources_with_det = len(all_detections)
            total_dets = sum(len(d) for d in all_detections.values())
            print(f"\n[{time.time()-t0:.1f}s] Initial scan complete: {n_sources_with_det} sources with "
                  f"{total_dets} total line detections", flush=True)

            # Tally and save to cache
            line_counts = tally_detections(all_detections=all_detections)
            save_detections_cache(all_detections=all_detections, line_counts=line_counts)

        print("\nDetection counts per line (before verification):")
        for ln, cnt in sorted(line_counts.items(), key=lambda x: -x[1]):
            print(f"  {ln}: {cnt}")

        # Verify detections by comparing to background
        # This step is optional but recommended - it can be slow
        if args.skip_verification:
            print(f"\n[{time.time()-t0:.1f}s] Skipping background verification step (--skip-verification)", flush=True)
            verified_detections = all_detections
        else:
            # Check for cached verified detections
            cached_verified = load_verified_cache()
            if cached_verified is not None and is_cache_valid():
                print(f"\n[{time.time()-t0:.1f}s] Using cached verified detections", flush=True)
                verified_detections = cached_verified
            else:
                print(f"\n[{time.time()-t0:.1f}s] Verifying detections against background annuli...", flush=True)
                print(f"[{time.time()-t0:.1f}s] (This step may extract background spectra from parent cubes - may be slow)", flush=True)
                verified_detections = verify_source_detections(
                    all_detections, catalog, linelist,
                    spectrum_dir=SPECTRUM_DIR, prefix=CATALOG_NAME_PREFIX
                )
                save_verified_cache(verified_detections)

                # Report verification statistics
                n_verified = len(verified_detections)
                total_verified = sum(len(d) for d in verified_detections.values())
                n_rejected_sources = len(all_detections) - n_verified
                total_rejected = total_dets - total_verified
                print(f"\n[{time.time()-t0:.1f}s] Verification complete:")
                print(f"[{time.time()-t0:.1f}s]   Verified: {n_verified} sources with {total_verified} detections")
                print(f"[{time.time()-t0:.1f}s]   Rejected: {n_rejected_sources} sources, {total_rejected} detections", flush=True)

        # Use verified detections for catalog augmentation and plotting
        final_detections = verified_detections
        final_line_counts = tally_detections(final_detections)
        
        # Keep all_detections for plotting (to show background-rejected lines)
        all_detections_for_plotting = all_detections

        print("\nFinal detection counts per line (after verification):")
        for ln, cnt in sorted(final_line_counts.items(), key=lambda x: -x[1]):
            print(f"  {ln}: {cnt}")

    # Top sources by number of different lines detected
    n_lines_per_source = {src_id: len(dets)
                         for src_id, dets in final_detections.items()}
    top_sources = sorted(n_lines_per_source.items(),
                        key=lambda x: -x[1])[:40]

    print("\n" + "="*70)
    print("TOP 40 SOURCES BY NUMBER OF DIFFERENT LINES DETECTED")
    print("="*70)
    for rank, (src_id, n_lines) in enumerate(top_sources, 1):
        line_names = list(final_detections[src_id].keys())
        print(f"{rank:3d}. Source {src_id:5d}: {n_lines:2d} lines - "
              f"{', '.join(line_names[:5])}"
              f"{' ...' if len(line_names) > 5 else ''}")
    print("="*70)

    # Augment catalog
    print(f"\n[{time.time()-t0:.1f}s] Augmenting catalog...", flush=True)
    catalog, sparse_dets = augment_catalog(catalog=catalog, all_detections=final_detections,
                                           line_counts=final_line_counts)

    # Write outputs
    print(f"\n[{time.time()-t0:.1f}s] Writing augmented catalog to {OUTPUT_CATALOG_PATH}", flush=True)
    catalog.write(OUTPUT_CATALOG_PATH, overwrite=True)
    print(f"[{time.time()-t0:.1f}s] Writing augmented catalog to {OUTPUT_CATALOG_ECSV_PATH}",
          flush=True)
    catalog.write(OUTPUT_CATALOG_ECSV_PATH, format='ascii.ecsv',
                  overwrite=True)

    if sparse_dets:
        print(f"[{time.time()-t0:.1f}s] Writing {len(sparse_dets)} rare-line source entries to "
              f"{SPARSE_DETECTIONS_JSON}", flush=True)
        with open(SPARSE_DETECTIONS_JSON, 'w') as fh:
            json.dump(sparse_dets, fh, indent=2)

    # Generate plots
    print(f"\n[{time.time()-t0:.1f}s] " + "="*70)
    print(f"[{time.time()-t0:.1f}s] GENERATING PLOTS")
    print(f"[{time.time()-t0:.1f}s] " + "="*70)

    # Plot spectra for top 40 sources
    print(f"\n[{time.time()-t0:.1f}s] Creating spectral plots for top 40 sources in {SPECTRA_PLOT_DIR}/...", flush=True)
    for rank, (src_id, n_lines) in enumerate(top_sources, 1):
        if rank % 10 == 0:
            print(f"[{time.time()-t0:.1f}s]   Plotted {rank}/40 sources...", flush=True)
        plot_source_spectra_with_fits(source_id=src_id, all_detections=all_detections_for_plotting, catalog=catalog,
                                     linelist=linelist, spectrum_dir=SPECTRUM_DIR,
                                     prefix=CATALOG_NAME_PREFIX, output_dir=SPECTRA_PLOT_DIR,
                                     verified_detections=final_detections)
    print(f"[{time.time()-t0:.1f}s]   Spectral plots complete: {len(top_sources)} plots saved", flush=True)

    # Plot spatial distribution for each line
    print(f"\n[{time.time()-t0:.1f}s] Creating spatial plots for each line in {SPATIAL_PLOT_DIR}/...", flush=True)
    unique_lines = sorted(final_line_counts.keys())
    for idx, line_name in enumerate(unique_lines, 1):
        if idx % 10 == 0:
            print(f"[{time.time()-t0:.1f}s]   Plotted {idx}/{len(unique_lines)} lines...", flush=True)
        plot_spatial_distribution_by_line(line_name=line_name, all_detections=final_detections,
                                         catalog=catalog, linelist=linelist, output_dir=SPATIAL_PLOT_DIR)
    print(f"[{time.time()-t0:.1f}s]   Spatial line plots complete: {len(unique_lines)} plots saved", flush=True)

    # Plot spatial distribution by number of lines
    print(f"\n[{time.time()-t0:.1f}s] Creating spatial plot by number of lines...", flush=True)
    plot_spatial_by_n_lines(all_detections=final_detections, catalog=catalog, output_dir=SPATIAL_PLOT_DIR)
    print(f"[{time.time()-t0:.1f}s]   Spatial plot by N lines saved", flush=True)

    print(f"\n[{time.time()-t0:.1f}s] Done. Total elapsed time: {time.time()-t0:.1f}s ({(time.time()-t0)/60:.1f}m)", flush=True)
