"""
Write a script that will measure the fluxes within the regions specified by aces/data/regions/ACES_point_source_selection_20250415.reg in both the mosaic fields and the individual MOUS images
"""

import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import mad_std
from regions import Regions, CirclePixelRegion, CircleSkyRegion
from pathlib import Path
from aces import conf
import glob
import warnings
from astropy.table import Table
from spectral_cube import SpectralCube, Slice
from astropy import units as u
from radio_beam import Beam
from spectral_cube.utils import NoBeamError
from astropy.coordinates import SkyCoord
from tqdm.auto import tqdm
import re
from astropy.io.registry import IORegistryError


def suppress_wcs_warnings():
    """
    Suppress specific WCS warnings from spectral_cube and astropy.wcs
    """
    warnings.filterwarnings('ignore', category=Warning, message='.*WCS.*')
    warnings.filterwarnings('ignore', message='.*missing card.*')
    warnings.filterwarnings('ignore', module='spectral_cube.wcs_utils')
    warnings.filterwarnings('ignore', module='astropy.wcs')


def suppress_numpy_warnings():
    """
    Suppress specific NumPy warnings related to empty slices and invalid values
    """
    warnings.filterwarnings('ignore', category=RuntimeWarning, message='Mean of empty slice')
    warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered in divide')


def read_file(filename):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        suppress_wcs_warnings()
        suppress_numpy_warnings()
        try:
            return SpectralCube.read(filename)[0]
        except IORegistryError:
            return SpectralCube.read(filename, format='casa_image')[0]
        except Exception as ex:
            return Slice.from_hdu(fits.open(filename))


def create_beam_aperture(point_region, beam):
    """
    Create a circular aperture around a point source with radius equal to the beam major axis
    
    Parameters
    ----------
    point_region : regions.Region
        The point source region
    beam : radio_beam.Beam
        The beam object
        
    Returns
    -------
    aperture : regions.CircleSkyRegion
        A circular region with radius equal to the beam major axis
    """
    # Get the center position from the point region
    center = point_region.center
    
    # Create a circular aperture with radius equal to the beam major axis
    radius = beam.major
    aperture = CircleSkyRegion(center=center, radius=radius)
    
    return aperture


def get_flux_in_region(fitsfile, region, rms_region=None):
    """
    Measure the flux within a region in a FITS image
    
    Parameters
    ----------
    fitsfile : str
        Path to the FITS file
    region : regions.Region
        Astropy region object (point source)
    rms_region : regions.Region, optional
        Region to measure the RMS noise
        
    Returns
    -------
    flux : float
        Total flux in the region in Jy
    peak : float
        Peak flux in the region in Jy/beam
    rms : float
        RMS noise level in Jy/beam
    """
    
    with warnings.catch_warnings():
        suppress_wcs_warnings()
        suppress_numpy_warnings()
        
        data = cube = read_file(fitsfile)

        beam = cube.beam
        aperture = create_beam_aperture(region, beam)

        header = cube.header
        wcs = cube.wcs.celestial
        
        # Convert sky aperture to pixel coordinates
        pixel_region = aperture.to_pixel(wcs)
        mask = pixel_region.to_mask()
        
        # Get the cutout of the region
        cutout_data = mask.multiply(data)
        if cutout_data is None:
            if wcs.footprint_contains(region.center):
                warnings.warn(f"Footprint contained region {region} but cutout was blank")
            return np.nan, np.nan, np.nan
        
        # Calculate flux (sum of pixels)
        pixel_scale = wcs.proj_plane_pixel_area()
        beam_area = beam.sr.to(u.steradian)
        pixels_per_beam = (beam_area / (pixel_scale * np.pi / (4 * np.log(2)))).decompose()
        
        # Get non-nan pixels within the aperture
        valid_data = cutout_data[np.isfinite(cutout_data)]
        if len(valid_data) == 0:
            return np.nan, np.nan, np.nan
        
        # Total flux is sum of pixels divided by pixels per beam
        flux = np.sum(valid_data) / pixels_per_beam
        peak = np.max(valid_data)
        
        # Calculate RMS from rms_region if provided, otherwise from the data edges
        if rms_region is not None:
            rms_pixel_region = rms_region.to_pixel(wcs)
            rms_mask = rms_pixel_region.to_mask()
            rms_data = rms_mask.multiply(data)
            valid_rms_data = rms_data[np.isfinite(rms_data)]
            if len(valid_rms_data) > 0:
                rms = mad_std(valid_rms_data)
            else:
                rms = np.nan
        else:
            # Use the edges of the image for RMS calculation
            edge_width = 20
            edges = np.concatenate([
                data[:edge_width, :].ravel(),
                data[-edge_width:, :].ravel(),
                data[:, :edge_width].ravel(),
                data[:, -edge_width:].ravel()
            ])
            valid_edges = edges[np.isfinite(edges)]
            if len(valid_edges) > 0:
                rms = mad_std(valid_edges)
            else:
                rms = np.nan
            
        return flux, peak, rms

def plot_source_fluxes(results_table, output_dir):
    """
    Plot the flux extracted from each source vs the image it is extracted from.
    
    Parameters
    ----------
    results_table : astropy.table.Table
        Table containing the flux measurement results
    output_dir : str
        Directory to save the plot files
    """
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    # Make sure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a figure for all sources
    all_sources_pdf = os.path.join(output_dir, 'all_point_source_fluxes.pdf')
    
    # Get unique sources
    sources = sorted(set(results_table['region']))
    
    # Get unique image types
    image_types = sorted(set(results_table['type']))
    
    # Filter the results table to include only finite, nonzero flux measurements
    valid_results = results_table[np.isfinite(results_table['flux_jy']) & (results_table['flux_jy'] > 0)]
    
    # If no valid results, return early
    if len(valid_results) == 0:
        print("No valid flux measurements found.")
        return
    
    # Get unique filenames for visualization, only those that have valid measurements
    unique_files = sorted(set(valid_results['filename']))
    
    # Strip suffixes from filenames for cleaner display
    # This will remove .fits, .image.pbcor.fits, etc.
    stripped_files = [re.sub(r'\.(?:fits|image\.tt0\.pbcor\.fits|image\.pbcor\.fits|tt0\.pbcor\.fits)$', '', fname) 
                     for fname in unique_files]
    
    file_to_idx = {fname: idx for idx, fname in enumerate(unique_files)}
    
    # Setup colors for sources
    colors = plt.cm.tab20(np.linspace(0, 1, len(sources)))
    source_to_color = {source: colors[i] for i, source in enumerate(sources)}
    
    # Setup markers for image types
    markers = {'mosaic': 'o', 'mous': 's'}
    
    with PdfPages(all_sources_pdf) as pdf:
        # 1. Plot fluxes for all sources across all files
        plt.figure(figsize=(14, 10))
        
        # Track which x positions have data for proper labeling
        used_x_positions = []
        
        for source in sources:
            source_mask = valid_results['region'] == source
            for img_type in image_types:
                type_mask = valid_results['type'] == img_type
                subset = valid_results[source_mask & type_mask]
                
                if len(subset) > 0:
                    # Convert filenames to numeric indices for x-axis
                    x_values = [file_to_idx[fname] for fname in subset['filename']]
                    used_x_positions.extend(x_values)
                    
                    plt.scatter(
                        x_values, 
                        subset['flux_jy'], 
                        label=f"{source} ({img_type})",
                        color=source_to_color[source],
                        marker=markers[img_type],
                        alpha=0.7,
                        edgecolor='k'
                    )
        
        # Only add labels for x positions that have data
        used_x_positions = sorted(set(used_x_positions))
        plt.xticks(used_x_positions, [stripped_files[unique_files.index(unique_files[i])] for i in used_x_positions], rotation=90)
        
        plt.xlabel('Image File')
        plt.ylabel('Flux (Jy)')
        plt.title('Flux Measurements for All Point Sources Across All Files')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        pdf.savefig()
        plt.close()
        
        # 2. Plot fluxes for each source individually
        for source in sources:
            source_mask = valid_results['region'] == source
            source_data = valid_results[source_mask]
            
            if len(source_data) == 0:
                continue
                
            plt.figure(figsize=(14, 10))
            used_x_positions = []
            
            for img_type in image_types:
                type_mask = source_data['type'] == img_type
                subset = source_data[type_mask]
                
                if len(subset) > 0:
                    # Convert filenames to numeric indices for x-axis
                    x_values = [file_to_idx[fname] for fname in subset['filename']]
                    used_x_positions.extend(x_values)
                    
                    plt.scatter(
                        x_values, 
                        subset['flux_jy'], 
                        label=f"{img_type}",
                        marker=markers[img_type],
                        alpha=0.7,
                        edgecolor='k'
                    )
                    
                    # Add errorbars based on RMS
                    plt.errorbar(
                        x_values, 
                        subset['flux_jy'], 
                        yerr=subset['rms_jy_beam'],
                        fmt='none',
                        alpha=0.5
                    )
            
            # Only add labels for x positions that have data
            used_x_positions = sorted(set(used_x_positions))
            plt.xticks(used_x_positions, [stripped_files[unique_files.index(unique_files[i])] for i in used_x_positions], rotation=90)
            
            # Get coordinates for the source
            if len(source_data) > 0:
                coords = source_data['coords'][0]  # Take the first entry's coordinates
            else:
                coords = "Unknown"
                
            plt.xlabel('Image File')
            plt.ylabel('Flux (Jy)')
            plt.title(f'Flux Measurements for Source {source} (Coordinates: {coords})')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.subplots_adjust(bottom=0.3)  # Make room for rotated labels
            pdf.savefig()
            plt.close()
            
        # 3. Plot a comparison of mosaic vs MOUS fluxes
        plt.figure(figsize=(10, 10))
        
        # For each source, compare average mosaic flux vs average MOUS flux
        mosaic_fluxes = []
        mous_fluxes = []
        source_labels = []
        
        for source in sources:
            source_mask = valid_results['region'] == source
            
            # Get average mosaic flux
            mosaic_mask = (valid_results['type'] == 'mosaic') & source_mask
            if np.sum(mosaic_mask) > 0:
                mosaic_flux = np.mean(valid_results['flux_jy'][mosaic_mask])
            else:
                mosaic_flux = np.nan
                
            # Get average MOUS flux
            mous_mask = (valid_results['type'] == 'mous') & source_mask
            if np.sum(mous_mask) > 0:
                mous_flux = np.mean(valid_results['flux_jy'][mous_mask])
            else:
                mous_flux = np.nan
                
            if not (np.isnan(mosaic_flux) or np.isnan(mous_flux)):
                mosaic_fluxes.append(mosaic_flux)
                mous_fluxes.append(mous_flux)
                source_labels.append(source)
        
        if len(mosaic_fluxes) > 0:
            # Plot the comparison
            scatter = plt.scatter(mous_fluxes, mosaic_fluxes, c=range(len(source_labels)), 
                      cmap='viridis', s=100, edgecolor='k')
            
            # Add perfect correlation line
            max_val = max(max(mosaic_fluxes), max(mous_fluxes))
            min_val = min(min(mosaic_fluxes), min(mous_fluxes))
            unity_line = plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='1:1 Line')
            
            # Add labels for each point
            for i, txt in enumerate(source_labels):
                plt.annotate(txt, (mous_fluxes[i], mosaic_fluxes[i]), 
                             xytext=(5, 5), textcoords='offset points')
            
            # Add colorbar as a legend for sources
            cbar = plt.colorbar(scatter, label='Source Index')
            cbar.set_ticks(range(len(source_labels)))
            cbar.set_ticklabels(source_labels)
            
            # Add legend for unity line
            plt.legend(handles=unity_line, loc='lower right')
            
            plt.xlabel('Average MOUS Flux (Jy)')
            plt.ylabel('Average Mosaic Flux (Jy)')
            plt.title('Comparison of Mosaic vs MOUS Fluxes')
            plt.grid(True, alpha=0.3)
            plt.axis('equal')
            plt.tight_layout()
            pdf.savefig()
            plt.close()
            
        # 4. NEW: Plot individual mosaic vs MOUS comparisons by spectral window groups
        # Define spectral window group patterns
        spw_groups = {
            'spw33_35': {
                'mosaic_pattern': r'spw33_35.*mosaic',
                'mous_pattern': r'spw33_35'
            },
            'spw25_27': {
                'mosaic_pattern': r'spw25_27.*mosaic',
                'mous_pattern': r'spw25_27'
            },
            'aggregate': {
                'mosaic_pattern': r'12m_continuum_commonbeam_circular_reimaged_mosaic',
                'mous_pattern': r'spw25_27_29_31_33_35'
            }
        }
        
        # Get unique mosaic files
        mosaic_files = sorted(set(valid_results[valid_results['type'] == 'mosaic']['filename']))
        
        # Get all MOUS files (excluding v1)
        all_mous_files = sorted(set(valid_results[valid_results['type'] == 'mous']['filename']))
        all_mous_files = [f for f in all_mous_files if 'v1' not in f]
        
        # For each spectral window group
        for group_name, patterns in spw_groups.items():
            print(f"Creating plots for {group_name} group")
            
            # Filter mosaic files for this group
            group_mosaic_files = [f for f in mosaic_files if re.search(patterns['mosaic_pattern'], f)]
            if not group_mosaic_files:
                print(f"No mosaic files found matching pattern {patterns['mosaic_pattern']}")
                continue
                
            # Filter MOUS files for this group
            group_mous_files = [f for f in all_mous_files if re.search(patterns['mous_pattern'], f)]
            if not group_mous_files:
                print(f"No MOUS files found matching pattern {patterns['mous_pattern']}")
                continue
                
            print(f"Found {len(group_mosaic_files)} mosaic files and {len(group_mous_files)} MOUS files for {group_name}")
            
            # For each mosaic file in this group
            for mosaic_file in group_mosaic_files:
                # Strip the suffix for display
                stripped_mosaic = re.sub(r'\.(?:fits|image\.tt0\.pbcor\.fits|image\.pbcor\.fits|tt0\.pbcor\.fits)$', '', mosaic_file)
                
                # Create a figure with two subplots side by side
                fig, (ax_flux, ax_peak) = plt.subplots(1, 2, figsize=(20, 10))
                
                # Data for the plot
                mosaic_data = valid_results[(valid_results['filename'] == mosaic_file) & 
                                           (valid_results['type'] == 'mosaic')]
                
                # Set up colors for MOUS files
                mous_colors = plt.cm.tab20(np.linspace(0, 1, len(group_mous_files)))
                mous_to_color = {fname: mous_colors[i] for i, fname in enumerate(group_mous_files)}
                
                # Keep track of which MOUS files are actually plotted
                plotted_mous_files = set()
                
                # Collect data grouped by source
                for source in sources:
                    mosaic_source_mask = (mosaic_data['region'] == source)
                    
                    # Skip if no mosaic data for this source
                    if np.sum(mosaic_source_mask) == 0:
                        continue
                        
                    mosaic_flux = mosaic_data['flux_jy'][mosaic_source_mask][0]  # Single value for this mosaic
                    mosaic_peak = mosaic_data['peak_jy_beam'][mosaic_source_mask][0]
                    
                    # Skip if mosaic flux or peak is zero or negative
                    if mosaic_flux <= 0 or mosaic_peak <= 0:
                        continue
                    
                    # Get all MOUS data for this source, but only from the current spectral window group
                    mous_source_data = valid_results[(valid_results['type'] == 'mous') & 
                                                   (valid_results['region'] == source)]
                    
                    # Filter to only include MOUS files in this group
                    mous_source_data = mous_source_data[[f in group_mous_files for f in mous_source_data['filename']]]
                    
                    # Skip if no MOUS data for this source
                    if len(mous_source_data) == 0:
                        continue
                    
                    # For each MOUS file with this source
                    for mous_file in set(mous_source_data['filename']):
                        # Skip if this MOUS file isn't in our filtered list
                        if mous_file not in group_mous_files:
                            continue
                            
                        mous_file_mask = mous_source_data['filename'] == mous_file
                        mous_flux = mous_source_data['flux_jy'][mous_file_mask][0]
                        mous_peak = mous_source_data['peak_jy_beam'][mous_file_mask][0]
                        
                        # Skip if MOUS flux or peak is zero or negative
                        if mous_flux <= 0 or mous_peak <= 0:
                            continue
                        
                        # Strip suffix for display
                        stripped_mous = re.sub(r'\.(?:fits|image\.tt0\.pbcor\.fits|image\.pbcor\.fits|tt0\.pbcor\.fits)$', '', mous_file)
                        
                        # Add to our list of plotted MOUS files
                        plotted_mous_files.add(mous_file)
                        
                        # Plot flux comparison
                        ax_flux.scatter(
                            mous_flux,
                            mosaic_flux,
                            color=mous_to_color[mous_file],
                            alpha=0.7,
                            s=50,
                            label=stripped_mous,
                            marker='o'
                        )
                        
                        # Annotate with source name
                        ax_flux.annotate(
                            source,
                            (mous_flux, mosaic_flux),
                            xytext=(5, 5),
                            textcoords='offset points',
                            fontsize=8
                        )
                        
                        # Plot peak comparison
                        ax_peak.scatter(
                            mous_peak,
                            mosaic_peak,
                            color=mous_to_color[mous_file],
                            alpha=0.7,
                            s=50,
                            label=stripped_mous,
                            marker='o'
                        )
                        
                        # Annotate with source name
                        ax_peak.annotate(
                            source,
                            (mous_peak, mosaic_peak),
                            xytext=(5, 5),
                            textcoords='offset points',
                            fontsize=8
                        )
                
                # If no MOUS files were plotted, skip this mosaic
                if not plotted_mous_files:
                    plt.close(fig)
                    continue
                    
                # Add 1:1 line to flux plot
                x_range = ax_flux.get_xlim()
                y_range = ax_flux.get_ylim()
                min_val = min(x_range[0], y_range[0])
                max_val = max(x_range[1], y_range[1])
                unity_line = ax_flux.plot([min_val, max_val], [min_val, max_val], 'k--', 
                                     alpha=0.5, label='1:1 Line')
                
                # Add 1:1 line to peak plot
                x_range = ax_peak.get_xlim()
                y_range = ax_peak.get_ylim()
                min_val = min(x_range[0], y_range[0])
                max_val = max(x_range[1], y_range[1])
                ax_peak.plot([min_val, max_val], [min_val, max_val], 'k--', 
                            alpha=0.5, label='1:1 Line')
                
                # Configure the flux subplot
                ax_flux.set_xlabel('MOUS Integrated Flux (Jy)')
                ax_flux.set_ylabel(f'Mosaic Integrated Flux from {stripped_mosaic} (Jy)')
                ax_flux.set_title(f'Integrated Flux Comparison')
                ax_flux.grid(True, alpha=0.3)
                
                # Configure the peak subplot
                ax_peak.set_xlabel('MOUS Peak Intensity (Jy/beam)')
                ax_peak.set_ylabel(f'Mosaic Peak Intensity from {stripped_mosaic} (Jy/beam)')
                ax_peak.set_title(f'Peak Intensity Comparison')
                ax_peak.grid(True, alpha=0.3)
                
                # Create a legend with the MOUS files that were actually plotted
                # Get unique labels to avoid duplicates in legend
                handles, labels = [], []
                for mous_file in sorted(plotted_mous_files):
                    stripped_mous = re.sub(r'\.(?:fits|image\.tt0\.pbcor\.fits|image\.pbcor\.fits|tt0\.pbcor\.fits)$', '', mous_file)
                    handle = plt.Line2D([], [], marker='o', color=mous_to_color[mous_file], 
                                       markersize=8, linestyle='None', label=stripped_mous)
                    handles.append(handle)
                    labels.append(stripped_mous)
                
                # Add unity line to legend
                unity_handle = plt.Line2D([], [], linestyle='--', color='k', alpha=0.5, label='1:1 Line')
                handles.append(unity_handle)
                labels.append('1:1 Line')
                
                # Add the legend to the figure
                # If we have many MOUS files, make the legend smaller and use more columns
                n_cols = min(5, max(3, len(handles) // 3))
                fontsize = 'medium' if len(handles) < 15 else 'small'
                
                # Position the legend higher on the page
                fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.15),
                          ncol=n_cols, frameon=True, fontsize=fontsize)
                
                # Add an overall title including the group name
                fig.suptitle(f'{group_name.upper()}: Comparison of {stripped_mosaic} vs Individual MOUS Measurements', fontsize=16)
                
                plt.tight_layout()
                # Adjust bottom margin to make room for the legend
                fig.subplots_adjust(bottom=0.3)  # Make more room for the legend
                pdf.savefig(fig)
                plt.close(fig)
    
    print(f"Plots saved to {all_sources_pdf}")
    
    # Also save a CSV file with a summary of the average fluxes per source
    summary_data = []
    for source in sources:
        source_mask = valid_results['region'] == source
        
        # Get average mosaic flux
        mosaic_mask = (valid_results['type'] == 'mosaic') & source_mask
        if np.sum(mosaic_mask) > 0:
            mosaic_flux = np.mean(valid_results['flux_jy'][mosaic_mask])
            mosaic_std = np.std(valid_results['flux_jy'][mosaic_mask])
        else:
            mosaic_flux = np.nan
            mosaic_std = np.nan
            
        # Get average MOUS flux
        mous_mask = (valid_results['type'] == 'mous') & source_mask
        if np.sum(mous_mask) > 0:
            mous_flux = np.mean(valid_results['flux_jy'][mous_mask])
            mous_std = np.std(valid_results['flux_jy'][mous_mask])
        else:
            mous_flux = np.nan
            mous_std = np.nan
            
        # Get coordinates for the source
        if np.sum(source_mask) > 0:
            coords = valid_results['coords'][np.where(source_mask)[0][0]]
        else:
            coords = "Unknown"
            
        summary_data.append({
            'source': source,
            'coords': coords,
            'mosaic_flux_jy': mosaic_flux,
            'mosaic_flux_std': mosaic_std,
            'mous_flux_jy': mous_flux,
            'mous_flux_std': mous_std,
            'flux_ratio': mosaic_flux / mous_flux if not (np.isnan(mosaic_flux) or np.isnan(mous_flux)) else np.nan
        })
    
    summary_table = Table(summary_data)
    summary_file = os.path.join(output_dir, 'point_source_flux_summary.ecsv')
    summary_table.write(summary_file, format='ascii.ecsv', overwrite=True)
    print(f"Summary table saved to {summary_file}")


def main():
    # Suppress WCS and NumPy warnings globally at start of main
    suppress_wcs_warnings()
    suppress_numpy_warnings()
    
    # Define paths
    basepath = conf.basepath
    region_file = os.path.join(basepath, 'reduction_ACES/aces/data/regions/ACES_point_source_selection_20250415.reg')
    
    # Check if region file exists
    if not os.path.exists(region_file):
        raise FileNotFoundError(f"Region file {region_file} not found!")
    
    # Load regions
    regions = Regions.read(region_file)
    
    # Find mosaic and MOUS images
    mosaic_pattern = os.path.join(basepath, 'mosaics', 'continuum', '*mosaic.fits')
    
    mosaic_files = glob.glob(mosaic_pattern)
    mosaic_files = [x for x in mosaic_files if ('tt1' not in x) and ('alpha' not in x) and ('rms' not in x) and ('weight' not in x)]
    assert len(mosaic_files) > 0, f"No mosaic files found in {mosaic_pattern}"
    
    # Initialize results table
    results = []
    
    print("Starting mosaic loop", flush=True)
    # Process mosaic files
    for fitsfile in tqdm(mosaic_files):
        try:
            beam = read_file(fitsfile).beam
            print(f"Processing mosaic: {os.path.basename(fitsfile)}", flush=True)
            print(f"Beam size: {beam.major.to(u.arcsec):.2f} x {beam.minor.to(u.arcsec):.2f}")
        except NoBeamError as e:
            # we can't measure fluxes from non-convolved ones, so we skip them
            continue
        except Exception as e:
            raise e
            # warnings.warn(f"Could not get beam from {fitsfile}: {str(e)}")
            continue
            
        for reg in regions:
            flux, peak, rms = get_flux_in_region(fitsfile, reg)
            results.append({
                'filename': os.path.basename(fitsfile),
                'type': 'mosaic',
                'region': reg.meta.get('text', 'unnamed'),
                'coords': f"{reg.center.galactic.l.deg:.6f},{reg.center.galactic.b.deg:.6f}",
                'flux_jy': flux,
                'peak_jy_beam': peak,
                'rms_jy_beam': rms,
                'beam_maj_arcsec': beam.major.to(u.arcsec).value,
                'beam_min_arcsec': beam.minor.to(u.arcsec).value,
                'beam_pa_deg': beam.pa.to(u.deg).value
            })
    assert len(results) > 0
    
    mous_pattern = os.path.join(basepath, 'data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001*/calibrated/working/*image.tt0.pbcor')  # Adjust pattern based on MOUS file location
    mous_files = glob.glob(mous_pattern)
    assert len(mous_files) > 40

    print(f"Starting MOUS loop with {len(mous_files)}", flush=True)
    # Process MOUS files
    for imagefile in tqdm(mous_files):
        print(f"Processing MOUS: {os.path.basename(imagefile)}")
        try:
            beam = read_file(imagefile).beam
            print(f"Beam size: {beam.major.to(u.arcsec):.2f} x {beam.minor.to(u.arcsec):.2f}")
        except Exception as e:
            warnings.warn(f"Could not get beam from {imagefile}: {str(e)}")
            continue
            
        for reg in regions:
            flux, peak, rms = get_flux_in_region(imagefile, reg)
            try:
                results.append({
                    'filename': os.path.basename(imagefile),
                    'type': 'mous',
                    'region': reg.meta.get('text', 'unnamed'),
                    'coords': f"{reg.center.galactic.l.deg:.6f},{reg.center.galactic.b.deg:.6f}",
                    'flux_jy': flux,
                    'peak_jy_beam': peak,
                    'rms_jy_beam': rms,
                    'beam_maj_arcsec': beam.major.to(u.arcsec).value,
                    'beam_min_arcsec': beam.minor.to(u.arcsec).value,
                    'beam_pa_deg': beam.pa.to(u.deg).value
                })
            except Exception as e:
                warnings.warn(f"Error processing region in {fitsfile}: {str(e)}")
    
    # Convert results to table and save
    results_table = Table(results)
    
    # Convert any Quantity objects to plain values before writing
    for column in results_table.colnames:
        if any(isinstance(val, u.Quantity) for val in results_table[column] if val is not None):
            results_table[column] = [val.value if isinstance(val, u.Quantity) else val 
                                    for val in results_table[column]]
    
    output_file = os.path.join(basepath, 'tables', 'point_source_fluxes.ecsv')
    results_table.write(output_file, format='ascii.ecsv', overwrite=True)
    print(f"Results saved to {output_file}")
    
    # Generate plots of source fluxes
    plot_dir = os.path.join(basepath, 'figures', 'point_source_analysis')
    plot_source_fluxes(results_table, plot_dir)

if __name__ == '__main__':
    main()
