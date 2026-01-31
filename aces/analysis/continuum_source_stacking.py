"""
The "final" catalog of sources in ACES is here:
/orange/adamginsburg/ACES/upload/compact_cont_source_catalog/official_cont_catalog_files/ACES_catalog_v0_README.txt
/orange/adamginsburg/ACES/upload/compact_cont_source_catalog/official_cont_catalog_files/aces_compact_catalog_v0.fits

The Herschel column density map, at 36" resolution, is here:
/orange/adamginsburg/cmz/herschel/column_properunits_conv36_source_only.fits

the better map from Yping Tang is here:
/orange/adamginsburg/cmz/tangyping_1mm/product_pixel_14_STMB/p1_med.fits
and the temperature:
/orange/adamginsburg/cmz/tangyping_1mm/product_pixel_14_STMB/p2_med.fits

The cmzoom catalog is at:
/orange/adamginsburg/cmz/cmzoom/catalog/catalog_robust.fits

There are spectral index maps:
/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_manual_alpha_mosaic.fits
/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_manual_alpha_rms_mosaic.fits
/orange/adamginsburg/ACES/mosaics/continuum/alpha_aces_meerkat.fits
/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_alpha_mosaic.fits
/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_alphaerror_mosaic.fits

Using this 'final' catalog, first:
1. compute the local column density and temperature derived from the Herschel data and the Tang 1mm data.  Just use the central position and store that as a new column in a new table:
/orange/adamginsburg/ACES/tables/aces_compact_catalog_v0_withExtras.fits
2. Compute the alpha values and uncertainties on those alpha values from each of the alpha maps listed above, again just at the central position, and store those as new columns in the same table.  Also compute spectral index by adding in the CMZoom flux and measuring it using the 1mm/3mm line ratio.
3. Bin the sources by column density with a range 5e21-5e23 and compute mean properties: flux, R_pc, Mean_CS, Mean_MS, alpha (each version of alpha).  Weight the mean by signal-to-noise of the flux.
4. Stack the sources in each column density bin to produce average images.  i.e., at each position, cut out a small image (e.g., 50x50 pixels) around each source in that bin, and average them together to produce a stacked image for that bin.  Weight the stack by the inverse variance of the image at each position.
The image and uncertainty to use are:
/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits
/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_maskedrms_mosaic.fits


Then,
1. make a version of the CMZoom complete catalog (catalog_complete.fits instead of catalog_robust.fits) in which the flux and flux uncertainty at the position of the CMZoom source in the ACES data is recorded.  The ACES image and uncertainty are:
/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_maskedrms_mosaic.fits
/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits
2. Compute the spectral index between CMZoom and ACES and limits on the spectral index if there is a nondetection in ACES.
3. Create a stacked image of ACES at the positions of CMZoom sources with ACES nondetections.  Again, weight the stack by the inverse variance of the ACES image.


"""

#!/usr/bin/env python3
"""
ACES Catalog Enhancement and Stacking Analysis

This script performs comprehensive analysis on the ACES compact continuum source catalog:
1. Extract column density and temperature from Herschel and Tang 1mm data
2. Extract spectral index values from various alpha maps
3. Compute spectral indices from CMZoom flux (1mm/3mm ratio)
4. Bin sources by column density and compute weighted mean properties
5. Stack images in each column density bin
6. Analyze CMZoom complete catalog with ACES flux measurements
7. Compute spectral indices between CMZoom and ACES
8. Stack ACES images at CMZoom nondetection positions

Author: GitHub Copilot
Date: January 29, 2026
"""

from astropy.io import fits
import numpy as np
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import astropy.units as u
from scipy.stats import binned_statistic
from astropy.nddata.utils import NoOverlapError
import matplotlib.pyplot as plt
import warnings
import matplotlib.gridspec as gridspec
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ==================== FILE PATHS ====================
# ACES catalog
ACES_CATALOG_PATH = "/orange/adamginsburg/ACES/upload/compact_cont_source_catalog/official_cont_catalog_files/aces_compact_catalog_v0.fits"
OUTPUT_CATALOG_PATH = "/orange/adamginsburg/ACES/tables/aces_compact_catalog_v0_withExtras.fits"

# Herschel and Tang data
HERSCHEL_COLDENS_PATH = "/orange/adamginsburg/cmz/herschel/column_properunits_conv36_source_only.fits"
TANG_COLDENS_PATH = "/orange/adamginsburg/cmz/tangyping_1mm/product_pixel_14_STMB/p1_med.fits"
TANG_TEMP_PATH = "/orange/adamginsburg/cmz/tangyping_1mm/product_pixel_14_STMB/p2_med.fits"

# Spectral index maps
ALPHA_MAPS = {
    'alpha_manual': "/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_manual_alpha_mosaic.fits",
    'alpha_manual_rms': "/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_manual_alpha_rms_mosaic.fits",
    'alpha_aces_meerkat': "/orange/adamginsburg/ACES/mosaics/continuum/alpha_aces_meerkat.fits",
    'alpha_auto': "/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_alpha_mosaic.fits",
    'alpha_auto_error': "/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_alphaerror_mosaic.fits"
}

# ACES mosaic images
ACES_IMAGE_PATH = "/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic.fits"
ACES_RMS_PATH = "/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_maskedrms_mosaic.fits"

# CMZoom catalogs
CMZOOM_ROBUST_PATH = "/orange/adamginsburg/cmz/cmzoom/catalog/catalog_robust.fits"
CMZOOM_COMPLETE_PATH = "/orange/adamginsburg/cmz/cmzoom/catalog/catalog_complete.fits"


# ==================== UTILITY FUNCTIONS ====================

def extract_value_at_position(image_path, lon, lat, wcs=None):
    """
    Extract pixel value from a FITS image at given sky coordinates.

    Parameters
    ----------
    image_path : str
        Path to FITS image
    lon : float or array
        Galactic longitude(s) in degrees
    lat : float or array
        Galactic latitude(s) in degrees
    wcs : WCS object, optional
        Pre-loaded WCS. If None, will be loaded from the image.

    Returns
    -------
    values : float or array
        Extracted pixel value(s)
    """
    with fits.open(image_path) as hdul:
        data = hdul[0].data
        if wcs is None:
            wcs = WCS(hdul[0].header, naxis=2)

        # Handle different dimensionality
        if data.ndim > 2:
            data = data.squeeze()

        # Ensure WCS is 2D
        if wcs.naxis > 2:
            wcs = wcs.celestial

        # Convert coordinates to pixel positions
        # Handle case where lon/lat might already have units
        if hasattr(lon, 'unit') and lon.unit is not None:
            coords = SkyCoord(lon, lat, frame='galactic')
        else:
            coords = SkyCoord(lon, lat, unit='deg', frame='galactic')
        px, py = wcs.world_to_pixel(coords)

        # Handle scalar vs array
        if np.isscalar(px):
            px, py = int(round(px)), int(round(py))
            if 0 <= px < data.shape[1] and 0 <= py < data.shape[0]:
                return data[py, px]
            else:
                return np.nan
        else:
            values = np.full_like(px, np.nan)
            valid = (px >= 0) & (px < data.shape[1]) & (py >= 0) & (py < data.shape[0])
            px_int = np.round(px[valid]).astype(int)
            py_int = np.round(py[valid]).astype(int)
            values[valid] = data[py_int, px_int]
            return values


def compute_spectral_index(flux1, flux2, freq1, freq2, flux1_err=None, flux2_err=None):
    """
    Compute spectral index: S_nu ~ nu^alpha

    alpha = ln(S1/S2) / ln(nu1/nu2)

    Parameters
    ----------
    flux1, flux2 : float or array
        Flux densities at two frequencies
    freq1, freq2 : float
        Frequencies in GHz
    flux1_err, flux2_err : float or array, optional
        Uncertainties in flux densities

    Returns
    -------
    alpha : float or array
        Spectral index
    alpha_err : float or array (if errors provided)
        Uncertainty in spectral index
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        alpha = np.log(flux1 / flux2) / np.log(freq1 / freq2)

    if flux1_err is not None and flux2_err is not None:
        # Error propagation for spectral index
        with np.errstate(divide='ignore', invalid='ignore'):
            term1 = (flux1_err / flux1) ** 2
            term2 = (flux2_err / flux2) ** 2
            alpha_err = np.abs(alpha) * np.sqrt(term1 + term2)
        return alpha, alpha_err

    return alpha


def weighted_mean(values, weights):
    """
    Compute weighted mean, handling NaNs.

    Parameters
    ----------
    values : array
        Values to average
    weights : array
        Weights for each value

    Returns
    -------
    mean : float
        Weighted mean
    """
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if np.sum(mask) == 0:
        return np.nan
    return np.sum(values[mask] * weights[mask]) / np.sum(weights[mask])


# ==================== MAIN PROCESSING FUNCTIONS ====================

def add_column_density_and_temperature(catalog):
    """
    Add column density and temperature columns from Herschel and Tang data.

    Parameters
    ----------
    catalog : Table
        ACES catalog table

    Returns
    -------
    catalog : Table
        Updated catalog with new columns
    """
    print("Extracting column density and temperature...")

    lon = catalog['GLON_peak']
    lat = catalog['GLAT_peak']

    # Extract Herschel column density
    print("  - Herschel column density...")
    herschel_coldens = extract_value_at_position(HERSCHEL_COLDENS_PATH, lon, lat)
    catalog['herschel_coldens'] = Column(herschel_coldens, unit='cm-2',
                                         description='Herschel column density at source position')

    # Extract Tang column density
    print("  - Tang 1mm column density...")
    tang_coldens = extract_value_at_position(TANG_COLDENS_PATH, lon, lat)
    catalog['tang_coldens'] = Column(tang_coldens, unit='cm-2',
                                     description='Tang 1mm column density at source position')

    # Extract Tang temperature
    print("  - Tang temperature...")
    tang_temp = extract_value_at_position(TANG_TEMP_PATH, lon, lat)
    catalog['tang_temp'] = Column(tang_temp, unit='K',
                                  description='Tang 1mm temperature at source position')

    return catalog


def add_spectral_indices(catalog):
    """
    Add spectral index values from various alpha maps.

    Parameters
    ----------
    catalog : Table
        ACES catalog table

    Returns
    -------
    catalog : Table
        Updated catalog with spectral index columns
    """
    print("Extracting spectral indices from alpha maps...")

    lon = catalog['GLON_peak']
    lat = catalog['GLAT_peak']

    for key, path in ALPHA_MAPS.items():
        print(f"  - {key}...")
        alpha_values = extract_value_at_position(path, lon, lat)
        catalog[key] = Column(alpha_values, description=f'Spectral index from {key}')

    return catalog


def add_cmzoom_spectral_index(catalog):
    """
    Compute spectral index using CMZoom flux and ACES flux (1mm/3mm ratio).

    Parameters
    ----------
    catalog : Table
        ACES catalog table with cmzoom_id column

    Returns
    -------
    catalog : Table
        Updated catalog with CMZoom spectral index
    """
    print("Computing CMZoom-ACES spectral indices...")

    # Load CMZoom catalog
    cmzoom_cat = Table.read(CMZOOM_ROBUST_PATH)

    # Initialize columns
    n_sources = len(catalog)
    cmzoom_flux = np.full(n_sources, np.nan)
    cmzoom_flux_err = np.full(n_sources, np.nan)
    alpha_cmzoom_aces = np.full(n_sources, np.nan)
    alpha_cmzoom_aces_err = np.full(n_sources, np.nan)

    # Match ACES sources with CMZoom
    for i, cmzoom_id in enumerate(catalog['cmzoom_id']):
        if cmzoom_id >= 0:  # Valid CMZoom match
            cmzoom_match = cmzoom_cat[cmzoom_cat['index'] == cmzoom_id]
            if len(cmzoom_match) > 0:
                # Try common flux column names
                flux_col = None
                for col in ['flux_integrated', 'flux', 'total_flux', 'int_flux']:
                    if col in cmzoom_match.colnames:
                        flux_col = col
                        break

                if flux_col:
                    cmzoom_flux[i] = cmzoom_match[flux_col][0]
                    flux_err_col = f'{flux_col}_err'
                    if flux_err_col in cmzoom_match.colnames:
                        cmzoom_flux_err[i] = cmzoom_match[flux_err_col][0]

    # Compute spectral index (CMZoom at 1mm, ACES at 3mm)
    # alpha = log(S_3mm / S_1mm) / log(nu_3mm / nu_1mm)
    freq_aces = 96.0  # GHz (3mm)
    freq_cmzoom = 260.0  # GHz (1mm)

    mask = np.isfinite(cmzoom_flux) & (cmzoom_flux > 0)
    if np.sum(mask) > 0:
        alpha_cmzoom_aces[mask], alpha_cmzoom_aces_err[mask] = compute_spectral_index(
            catalog['flux'][mask],
            cmzoom_flux[mask],
            freq_aces,
            freq_cmzoom,
            catalog['flux_err'][mask],
            cmzoom_flux_err[mask]
        )

    catalog['cmzoom_flux'] = Column(cmzoom_flux, unit='Jy',
                                    description='CMZoom 1mm flux for matched sources')
    catalog['cmzoom_flux_err'] = Column(cmzoom_flux_err, unit='Jy',
                                        description='CMZoom 1mm flux uncertainty')
    catalog['alpha_cmzoom_aces'] = Column(alpha_cmzoom_aces,
                                          description='Spectral index from ACES (3mm) to CMZoom (1mm)')
    catalog['alpha_cmzoom_aces_err'] = Column(alpha_cmzoom_aces_err,
                                              description='Uncertainty in ACES-CMZoom spectral index')

    return catalog


def bin_and_compute_mean_properties(catalog, output_table_path=None):
    """
    Bin sources by column density and compute weighted mean properties.

    Parameters
    ----------
    catalog : Table
        ACES catalog with column density information
    output_table_path : str, optional
        Path to save binned statistics table

    Returns
    -------
    binned_stats : Table
        Table with mean properties for each column density bin
    """
    print("Binning sources by column density and computing mean properties...")

    # Define column density bins (5e21 to 5e23 cm^-2)
    coldens_bins = np.logspace(np.log10(5e21), np.log10(5e23), 11)  # 10 bins

    # Use Herschel column density for binning (more complete than Tang)
    coldens = catalog['herschel_coldens']

    # Compute signal-to-noise weights
    snr = catalog['flux'] / catalog['flux_err']
    weights = snr

    # Properties to compute means for
    properties = ['flux', 'R_pc', 'Mean_CS', 'Mean_MS',
                  'alpha_manual', 'alpha_aces_meerkat', 'alpha_auto', 'alpha_cmzoom_aces']

    # Initialize results
    bin_centers = []
    bin_counts = []
    mean_properties = {prop: [] for prop in properties}
    err_properties = {prop: [] for prop in properties}

    for i in range(len(coldens_bins) - 1):
        bin_mask = (coldens >= coldens_bins[i]) & (coldens < coldens_bins[i+1])
        n_in_bin = np.sum(bin_mask)

        bin_centers.append(np.sqrt(coldens_bins[i] * coldens_bins[i+1]))
        bin_counts.append(n_in_bin)

        if n_in_bin > 0:
            for prop in properties:
                if prop in catalog.colnames:
                    values = catalog[prop][bin_mask]
                    bin_weights = weights[bin_mask]
                    mean_val = weighted_mean(values, bin_weights)
                    mean_properties[prop].append(mean_val)

                    # Compute weighted standard error
                    valid = np.isfinite(values) & np.isfinite(bin_weights) & (bin_weights > 0)
                    if np.sum(valid) > 1:
                        weighted_var = np.sum(bin_weights[valid] * (values[valid] - mean_val)**2) / np.sum(bin_weights[valid])
                        err = np.sqrt(weighted_var / np.sum(valid))  # Standard error
                    else:
                        err = np.nan
                    err_properties[prop].append(err)
                else:
                    mean_properties[prop].append(np.nan)
                    err_properties[prop].append(np.nan)
        else:
            for prop in properties:
                mean_properties[prop].append(np.nan)
                err_properties[prop].append(np.nan)

    # Create output table
    binned_stats = Table()
    binned_stats['coldens_bin_center'] = Column(bin_centers, unit='cm-2')
    binned_stats['n_sources'] = Column(bin_counts)

    for prop in properties:
        if prop in catalog.colnames:
            unit = catalog[prop].unit if hasattr(catalog[prop], 'unit') else None
            binned_stats[f'mean_{prop}'] = Column(mean_properties[prop], unit=unit)
            binned_stats[f'err_{prop}'] = Column(err_properties[prop], unit=unit)

    print(f"  Created {len(coldens_bins)-1} bins with sources in each")

    if output_table_path:
        binned_stats.write(output_table_path, overwrite=True)
        print(f"  Saved binned statistics to {output_table_path}")

    return binned_stats


def stack_images_by_coldens_bin(catalog, output_dir='/orange/adamginsburg/ACES/analysis/stacked_images'):
    """
    Stack images in each column density bin.

    Parameters
    ----------
    catalog : Table
        ACES catalog with column density information
    output_dir : str
        Directory to save stacked images

    Returns
    -------
    stacked_images : dict
        Dictionary with bin indices as keys and stacked image arrays as values
    """
    print("Stacking images by column density bin...")

    import os
    os.makedirs(output_dir, exist_ok=True)

    # Define column density bins
    coldens_bins = np.logspace(np.log10(5e21), np.log10(5e23), 11)
    coldens = catalog['herschel_coldens']

    # Load ACES image and RMS
    with fits.open(ACES_IMAGE_PATH) as hdul:
        aces_data = hdul[0].data.squeeze()
        aces_header = hdul[0].header
        aces_wcs = WCS(aces_header).celestial

    with fits.open(ACES_RMS_PATH) as hdul:
        aces_rms = hdul[0].data.squeeze()

    # Cutout size (50x50 pixels)
    cutout_size = (50, 50)

    stacked_images = {}

    for i in range(len(coldens_bins) - 1):
        print(f"  Bin {i+1}/{len(coldens_bins)-1}...")

        bin_mask = (coldens >= coldens_bins[i]) & (coldens < coldens_bins[i+1])
        n_in_bin = np.sum(bin_mask)

        if n_in_bin == 0:
            print(f"    No sources in bin")
            continue

        # Get sources in this bin
        bin_sources = catalog[bin_mask]

        # Initialize stacked arrays
        stacked_flux = np.zeros(cutout_size)
        stacked_weights = np.zeros(cutout_size)
        n_cutouts = 0

        for source in bin_sources:
            try:
                # Create SkyCoord for source position
                source_coord = SkyCoord(source['GLON_peak'],
                                        source['GLAT_peak'],
                                        unit='deg',
                                        frame='galactic')

                # Create cutout
                cutout = Cutout2D(aces_data, source_coord, cutout_size, wcs=aces_wcs)
                cutout_rms = Cutout2D(aces_rms, source_coord, cutout_size, wcs=aces_wcs)

                # Weight by inverse variance
                with np.errstate(divide='ignore', invalid='ignore'):
                    inv_var = 1.0 / cutout_rms.data**2
                    inv_var[~np.isfinite(inv_var)] = 0

                # Add to stack
                stacked_flux += cutout.data * inv_var
                stacked_weights += inv_var
                n_cutouts += 1

            except Exception as e:
                # Skip sources that can't be cut out (e.g., near edges)
                # but we need to know what specific exceptions cause edges, we don't want to accidentally skip any so we're going to run until we catch an exception, then if it is caused by something we understand (like an edge case), we can catch that exception and continue.
                raise e

        # Normalize by weights
        with np.errstate(divide='ignore', invalid='ignore'):
            stacked_image = stacked_flux / stacked_weights
            stacked_image[stacked_weights == 0] = np.nan

        stacked_images[i] = stacked_image

        print(f"    Stacked {n_cutouts} sources")

        # Save stacked image
        if n_cutouts > 0:
            output_path = os.path.join(output_dir, f'stacked_coldens_bin{i:02d}.fits')

            # Create a simple header for the stacked image
            header = fits.Header()
            header['NAXIS'] = 2
            header['NAXIS1'] = cutout_size[0]
            header['NAXIS2'] = cutout_size[1]
            header['CRPIX1'] = cutout_size[0] / 2
            header['CRPIX2'] = cutout_size[1] / 2
            header['NSOURCES'] = n_cutouts
            header['CDMIN'] = coldens_bins[i]
            header['CDMAX'] = coldens_bins[i+1]
            header['BUNIT'] = 'Jy/beam'
            header['COMMENT'] = f'Stacked image from {n_cutouts} sources'

            hdu = fits.PrimaryHDU(data=stacked_image, header=header)
            hdu.writeto(output_path, overwrite=True)
            print(f"    Saved to {output_path}")

    return stacked_images


def create_cmzoom_catalog_with_aces_flux(output_path='/orange/adamginsburg/ACES/tables/cmzoom_complete_with_aces.fits'):
    """
    Create CMZoom complete catalog with ACES flux measurements.

    Parameters
    ----------
    output_path : str
        Path to save the enhanced CMZoom catalog

    Returns
    -------
    catalog : Table
        CMZoom catalog with ACES flux measurements
    """
    print("Creating CMZoom catalog with ACES flux measurements...")

    # Load CMZoom complete catalog
    cmzoom_cat = Table.read(CMZOOM_COMPLETE_PATH)
    print(f"  Loaded {len(cmzoom_cat)} CMZoom sources")

    # Extract ACES flux and RMS at CMZoom positions
    # Assume CMZoom has GLON and GLAT columns (need to check actual column names)
    # Common column names might be: glon, glat, ra, dec, l, b, etc.

    # Check available columns
    print(f"  CMZoom catalog columns: {cmzoom_cat.colnames[:10]}...")  # Print first 10

    # Try to find coordinate columns
    lon_col = None
    lat_col = None
    for col in cmzoom_cat.colnames:
        col_lower = col.lower()
        if 'glon' in col_lower or col_lower in ['l', 'lon']:
            lon_col = col
        if 'glat' in col_lower or col_lower in ['b', 'lat']:
            lat_col = col

    if lon_col is None or lat_col is None:
        # Try RA/Dec instead
        for col in cmzoom_cat.colnames:
            col_lower = col.lower()
            if 'ra' in col_lower:
                lon_col = col
            if 'dec' in col_lower:
                lat_col = col

        if lon_col and lat_col:
            # Convert RA/Dec to Galactic
            coords = SkyCoord(cmzoom_cat[lon_col], cmzoom_cat[lat_col], unit='deg', frame='icrs')
            lon = coords.galactic.l.deg
            lat = coords.galactic.b.deg
        else:
            raise ValueError("Cannot find coordinate columns in CMZoom catalog")
    else:
        lon = cmzoom_cat[lon_col]
        lat = cmzoom_cat[lat_col]

    # Extract ACES flux and RMS
    print("  Extracting ACES flux at CMZoom positions...")
    aces_flux = extract_value_at_position(ACES_IMAGE_PATH, lon, lat)
    aces_rms = extract_value_at_position(ACES_RMS_PATH, lon, lat)

    # Add to catalog
    cmzoom_cat['aces_flux'] = Column(aces_flux, unit='Jy/beam',
                                     description='ACES 3mm flux at CMZoom position')
    cmzoom_cat['aces_rms'] = Column(aces_rms, unit='Jy/beam',
                                    description='ACES 3mm RMS at CMZoom position')

    # Compute spectral index between CMZoom and ACES
    print("  Computing CMZoom-ACES spectral indices...")

    # Get CMZoom flux
    cmzoom_flux_col = None
    for col in ['flux_integrated', 'flux', 'total_flux', 'peak_flux', 'int_flux']:
        if col in cmzoom_cat.colnames:
            cmzoom_flux_col = col
            break

    if cmzoom_flux_col is None:
        print("  Warning: Could not find flux column in CMZoom catalog")
        cmzoom_flux = np.full(len(cmzoom_cat), np.nan)
        cmzoom_flux_err = np.full(len(cmzoom_cat), np.nan)
    else:
        print(f"  Using flux column: {cmzoom_flux_col}")
        cmzoom_flux = cmzoom_cat[cmzoom_flux_col]
        flux_err_col = f'{cmzoom_flux_col}_err'
        if flux_err_col in cmzoom_cat.colnames:
            cmzoom_flux_err = cmzoom_cat[flux_err_col]
        else:
            # Estimate error from noise column if available
            if 'noise' in cmzoom_cat.colnames:
                cmzoom_flux_err = cmzoom_cat['noise']
            else:
                cmzoom_flux_err = np.full(len(cmzoom_cat), np.nan)

    # Frequencies
    freq_aces = 96.0  # GHz (3mm)
    freq_cmzoom = 226.0  # GHz (1mm)

    # Compute spectral index for detections
    alpha = np.full(len(cmzoom_cat), np.nan)
    alpha_err = np.full(len(cmzoom_cat), np.nan)
    alpha_limit = np.full(len(cmzoom_cat), np.nan)  # Upper/lower limits for nondetections

    # Detections: S/N > 5 in ACES
    detection_mask = (aces_flux / aces_rms) > 5

    if np.sum(detection_mask) > 0:
        alpha[detection_mask], alpha_err[detection_mask] = compute_spectral_index(
            aces_flux[detection_mask],
            cmzoom_flux[detection_mask],
            freq_aces,
            freq_cmzoom,
            aces_rms[detection_mask],
            cmzoom_flux_err[detection_mask]
        )

    # Non-detections: compute limit assuming 5-sigma upper limit
    nondetection_mask = ~detection_mask & np.isfinite(aces_rms) & (cmzoom_flux > 0)
    if np.sum(nondetection_mask) > 0:
        flux_limit = 5 * aces_rms[nondetection_mask]
        alpha_limit[nondetection_mask] = compute_spectral_index(
            flux_limit,
            cmzoom_flux[nondetection_mask],
            freq_aces,
            freq_cmzoom
        )

    cmzoom_cat['alpha_cmzoom_aces'] = Column(alpha,
                                             description='Spectral index from CMZoom to ACES')
    cmzoom_cat['alpha_cmzoom_aces_err'] = Column(alpha_err,
                                                 description='Uncertainty in spectral index')
    cmzoom_cat['alpha_limit'] = Column(alpha_limit,
                                       description='Spectral index limit for nondetections')
    cmzoom_cat['aces_detection'] = Column(detection_mask,
                                          description='ACES detection flag (S/N > 5)')

    # Save catalog
    cmzoom_cat.write(output_path, overwrite=True)
    print(f"  Saved enhanced CMZoom catalog to {output_path}")
    print(f"  ACES detections: {np.sum(detection_mask)}/{len(cmzoom_cat)}")

    return cmzoom_cat


def stack_aces_at_cmzoom_nondetections(cmzoom_catalog=None,
                                       cmzoom_catalog_path=None,
                                       output_path='/orange/adamginsburg/ACES/analysis/stacked_images/stacked_cmzoom_nondetections.fits',
                                       catalog_name='CMZoom'):
    """
    Stack ACES images at positions of CMZoom nondetections.

    Parameters
    ----------
    cmzoom_catalog : Table, optional
        CMZoom catalog with ACES measurements. If None, will load from file or create.
    cmzoom_catalog_path : str, optional
        Path to CMZoom catalog to load and add ACES flux to
    output_path : str
        Path to save stacked image
    catalog_name : str
        Name for display purposes

    Returns
    -------
    stacked_image : array
        Stacked ACES image
    """
    print(f"Stacking ACES images at {catalog_name} nondetection positions...")

    if cmzoom_catalog is None:
        if cmzoom_catalog_path is not None:
            # Load catalog and add ACES flux
            cmzoom_catalog = Table.read(cmzoom_catalog_path)
            print(f"  Loaded {len(cmzoom_catalog)} {catalog_name} sources")

            # Find coordinates
            lon_col = None
            lat_col = None
            for col in cmzoom_catalog.colnames:
                col_lower = col.lower()
                if 'glon' in col_lower or col_lower in ['l', 'lon']:
                    lon_col = col
                if 'glat' in col_lower or col_lower in ['b', 'lat']:
                    lat_col = col

            if lon_col is None or lat_col is None:
                for col in cmzoom_catalog.colnames:
                    col_lower = col.lower()
                    if 'ra' in col_lower:
                        lon_col = col
                    if 'dec' in col_lower:
                        lat_col = col

                if lon_col and lat_col:
                    coords = SkyCoord(cmzoom_catalog[lon_col], cmzoom_catalog[lat_col], unit='deg', frame='icrs')
                    lon = coords.galactic.l.deg
                    lat = coords.galactic.b.deg
            else:
                lon = cmzoom_catalog[lon_col]
                lat = cmzoom_catalog[lat_col]

            # Extract ACES flux and RMS
            print("  Extracting ACES flux at CMZoom positions...")
            aces_flux = extract_value_at_position(ACES_IMAGE_PATH, lon, lat)
            aces_rms = extract_value_at_position(ACES_RMS_PATH, lon, lat)

            # Add to catalog and compute detection flag
            cmzoom_catalog['aces_flux'] = Column(aces_flux, unit='Jy/beam')
            cmzoom_catalog['aces_rms'] = Column(aces_rms, unit='Jy/beam')
            detection_mask = (aces_flux / aces_rms) > 5
            cmzoom_catalog['aces_detection'] = Column(detection_mask)

            print(f"  ACES detections: {np.sum(detection_mask)}/{len(cmzoom_catalog)}")
        else:
            cmzoom_catalog = Table.read('/orange/adamginsburg/ACES/tables/cmzoom_complete_with_aces.fits')

    # Get nondetections
    nondetection_mask = cmzoom_catalog['aces_detection'] == False
    nondetections = cmzoom_catalog[nondetection_mask]

    print(f"  Found {len(nondetections)} nondetections to stack")

    if len(nondetections) == 0:
        print("  No nondetections to stack")
        return None

    # Load ACES image and RMS
    with fits.open(ACES_IMAGE_PATH) as hdul:
        aces_data = hdul[0].data.squeeze()
        aces_header = hdul[0].header
        aces_wcs = WCS(aces_header).celestial

    with fits.open(ACES_RMS_PATH) as hdul:
        aces_rms = hdul[0].data.squeeze()

    # Cutout size
    cutout_size = (50, 50)

    # Initialize stacked arrays
    stacked_flux = np.zeros(cutout_size)
    stacked_weights = np.zeros(cutout_size)
    n_cutouts = 0

    # Get coordinates
    lon_col, lat_col = None, None
    for col in cmzoom_catalog.colnames:
        col_lower = col.lower()
        if 'glon' in col_lower or col_lower in ['l', 'lon']:
            lon_col = col
        if 'glat' in col_lower or col_lower in ['b', 'lat']:
            lat_col = col

    for source in nondetections:
        try:
            if lon_col and lat_col:
                source_coord = SkyCoord(source[lon_col], source[lat_col], unit='deg', frame='galactic')
            else:
                # Try RA/Dec
                source_coord = SkyCoord(source['ra'], source['dec'], unit='deg', frame='icrs')

            # Create cutout
            cutout = Cutout2D(aces_data, source_coord, cutout_size, wcs=aces_wcs)
            cutout_rms = Cutout2D(aces_rms, source_coord, cutout_size, wcs=aces_wcs)

            # Weight by inverse variance
            with np.errstate(divide='ignore', invalid='ignore'):
                inv_var = 1.0 / cutout_rms.data**2
                inv_var[~np.isfinite(inv_var)] = 0

            # Add to stack
            to_stack = cutout.data * inv_var
            to_stack[~np.isfinite(to_stack)] = 0
            stacked_flux += to_stack
            stacked_weights += inv_var
            n_cutouts += 1

        except NoOverlapError:
            # skip because out of bounds - this case is OK.
            continue
        except ValueError as ex:
            if "operands could not be broadcast together" in str(ex):
                # skip because out of bounds - this case is OK.
                continue
            else:
                # re-raise
                raise ex
        except Exception as e:
            # as above, only catch exceptions we understand
            raise e

    # Normalize
    with np.errstate(divide='ignore', invalid='ignore'):
        stacked_image = stacked_flux / stacked_weights
        stacked_image[stacked_weights == 0] = np.nan

    print(f"  Stacked {n_cutouts} nondetections")

    # Save
    if n_cutouts > 0:
        import os
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        header = fits.Header()
        header['NAXIS'] = 2
        header['NAXIS1'] = cutout_size[0]
        header['NAXIS2'] = cutout_size[1]
        header['CRPIX1'] = cutout_size[0] / 2
        header['CRPIX2'] = cutout_size[1] / 2
        header['NSOURCES'] = n_cutouts
        header['BUNIT'] = 'Jy/beam'
        header['COMMENT'] = f'Stacked ACES at CMZoom nondetections'

        hdu = fits.PrimaryHDU(data=stacked_image, header=header)
        hdu.writeto(output_path, overwrite=True)
        print(f"  Saved to {output_path}")

    return stacked_image


def plot_binned_statistics(binned_stats, output_path='/orange/adamginsburg/ACES/diagnostic_plots/binned_statistics_plots.png'):
    """
    Create a grid of plots showing binned statistics vs column density.

    Parameters
    ----------
    binned_stats : Table
        Binned statistics table
    output_path : str
        Path to save the plot
    """
    print("\nCreating binned statistics plots...")

    # Define which properties to plot
    properties = ['n_sources', 'mean_flux', 'mean_R_pc', 'mean_Mean_CS', 'mean_Mean_MS',
                  'mean_alpha_manual', 'mean_alpha_aces_meerkat', 'mean_alpha_auto',
                  'mean_alpha_cmzoom_aces']

    # Filter to properties that exist and have at least some finite values
    properties_to_plot = []
    for prop in properties:
        if prop in binned_stats.colnames:
            if prop == 'n_sources' or np.sum(np.isfinite(binned_stats[prop])) > 0:
                properties_to_plot.append(prop)

    n_plots = len(properties_to_plot)

    # Create figure with grid layout
    ncols = 3
    nrows = int(np.ceil(n_plots / ncols))

    fig = plt.figure(figsize=(15, 4*nrows))
    gs = gridspec.GridSpec(nrows, ncols, figure=fig, hspace=0.3, wspace=0.3)

    coldens = binned_stats['coldens_bin_center']

    for idx, prop in enumerate(properties_to_plot):
        row = idx // ncols
        col = idx % ncols
        ax = fig.add_subplot(gs[row, col])

        values = binned_stats[prop]

        if prop == 'n_sources':
            # Histogram plot for number of sources
            ax.bar(np.arange(len(values)), values, width=0.8, alpha=0.7, edgecolor='black')
            ax.set_xlabel('Bin Number')
            ax.set_ylabel('Number of Sources')
            ax.set_title('Source Count by Bin')

            # Add total count
            total = np.sum(values)
            ax.text(0.95, 0.95, f'Total: {total}', transform=ax.transAxes,
                    ha='right', va='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        else:
            # Scatter plot with error bars for other properties
            mask = np.isfinite(values) & np.isfinite(coldens)

            # Get error bars if available
            err_prop = prop.replace('mean_', 'err_')
            if err_prop in binned_stats.colnames:
                errors = binned_stats[err_prop]
                err_mask = mask & np.isfinite(errors)
            else:
                errors = None
                err_mask = mask

            if np.sum(err_mask) > 0:
                if errors is not None:
                    ax.errorbar(coldens[err_mask], values[err_mask], yerr=errors[err_mask],
                                fmt='o-', markersize=8, linewidth=2, capsize=5, capthick=2, alpha=0.7)
                else:
                    ax.plot(coldens[mask], values[mask], 'o-', markersize=8, linewidth=2, alpha=0.7)

                ax.set_xscale('log')
                ax.set_xlabel('Column Density (cm$^{-2}$)')

                # Format y-label
                ylabel = prop.replace('mean_', '').replace('_', ' ').title()
                if 'alpha' in prop.lower():
                    ylabel = ylabel.replace('Alpha', 'α')
                ax.set_ylabel(ylabel)

                ax.set_title(ylabel)

                # Add horizontal line at y=0 for alpha plots
                if 'alpha' in prop.lower():
                    ax.axhline(0, color='k', linestyle='--', alpha=0.3, linewidth=1)
            else:
                ax.text(0.5, 0.5, 'No finite data', transform=ax.transAxes,
                        ha='center', va='center', fontsize=14)
                ax.set_xlabel('Column Density (cm$^{-2}$)')
                ax.set_ylabel(prop.replace('mean_', ''))

    # Add overall title
    fig.suptitle('ACES Binned Statistics vs Column Density', fontsize=16, y=0.995)

    # Save figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved plot to {output_path}")
    plt.close()

    return fig


# ==================== MAIN EXECUTION ====================

def main():
    """
    Main execution function.
    """
    print("="*70)
    print("ACES Catalog Enhancement and Stacking Analysis")
    print("="*70)

    # Load ACES catalog
    print("\nLoading ACES catalog...")
    catalog = Table.read(ACES_CATALOG_PATH)
    print(f"Loaded {len(catalog)} sources")

    # Part 1: Enhance ACES catalog
    print("\n" + "="*70)
    print("PART 1: Enhancing ACES Catalog")
    print("="*70)

    catalog = add_column_density_and_temperature(catalog)
    catalog = add_spectral_indices(catalog)
    catalog = add_cmzoom_spectral_index(catalog)

    # Save enhanced catalog
    print(f"\nSaving enhanced catalog to {OUTPUT_CATALOG_PATH}...")
    catalog.write(OUTPUT_CATALOG_PATH, overwrite=True)
    print("Done!")

    # Part 2: Bin and compute mean properties
    print("\n" + "="*70)
    print("PART 2: Binning by Column Density")
    print("="*70)

    binned_stats = bin_and_compute_mean_properties(
        catalog,
        output_table_path='/orange/adamginsburg/ACES/tables/aces_binned_by_coldens.fits'
    )
    print("\nBinned statistics:")
    print(binned_stats)

    # Plot binned statistics
    plot_binned_statistics(binned_stats)

    # Part 3: Stack images
    print("\n" + "="*70)
    print("PART 3: Stacking Images by Column Density")
    print("="*70)

    stacked_images = stack_images_by_coldens_bin(catalog)

    # Part 4: CMZoom analysis
    print("\n" + "="*70)
    print("PART 4: CMZoom Catalog Analysis")
    print("="*70)

    cmzoom_cat = create_cmzoom_catalog_with_aces_flux()

    # Part 5: Stack CMZoom nondetections
    print("\n" + "="*70)
    print("PART 5: Stacking CMZoom Nondetections")
    print("="*70)

    # Stack complete catalog nondetections (using the catalog we just created)
    print("\n--- CMZoom Complete Catalog ---")
    stack_aces_at_cmzoom_nondetections(
        cmzoom_catalog=cmzoom_cat,
        output_path='/orange/adamginsburg/ACES/analysis/stacked_images/stacked_cmzoom_complete_nondetections.fits',
        catalog_name='CMZoom Complete'
    )

    # Stack robust catalog nondetections (load from file)
    print("\n--- CMZoom Robust Catalog ---")
    stack_aces_at_cmzoom_nondetections(
        cmzoom_catalog_path=CMZOOM_ROBUST_PATH,
        output_path='/orange/adamginsburg/ACES/analysis/stacked_images/stacked_cmzoom_robust_nondetections.fits',
        catalog_name='CMZoom Robust'
    )

    print("\n" + "="*70)
    print("Analysis Complete!")
    print("="*70)


if __name__ == '__main__':
    main()
