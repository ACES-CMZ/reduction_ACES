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
from astropy.coordinates import SkyCoord


def suppress_wcs_warnings():
    """
    Suppress specific WCS warnings from spectral_cube and astropy.wcs
    """
    warnings.filterwarnings('ignore', category=Warning, message='.*WCS.*')
    warnings.filterwarnings('ignore', message='.*missing card.*')
    warnings.filterwarnings('ignore', module='spectral_cube.wcs_utils')
    warnings.filterwarnings('ignore', module='astropy.wcs')


def read_file(filename):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        suppress_wcs_warnings()
        try:
            return SpectralCube.read(filename)
        except Exception as ex:
            return Slice.from_hdu(filename)


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
    
    cube = read_file(fitsfile)

    if cube.ndim > 2:
        data = cube[0]
    else:
        data = cube

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
    beam_area = beam.sr.to(u.steradian).value
    pixels_per_beam = beam_area / (pixel_scale * np.pi / (4 * np.log(2)))
    
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
        rms = mad_std(rms_data[np.isfinite(rms_data)])
    else:
        # Use the edges of the image for RMS calculation
        edge_width = 20
        edges = np.concatenate([
            data[:edge_width, :].ravel(),
            data[-edge_width:, :].ravel(),
            data[:, :edge_width].ravel(),
            data[:, -edge_width:].ravel()
        ])
        rms = mad_std(edges[np.isfinite(edges)])
        
    return flux, peak, rms

def main():
    # Suppress WCS warnings globally at start of main
    suppress_wcs_warnings()
    
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
    
    # Initialize results table
    results = []
    
    print("Starting mosaic loop", flush=True)
    # Process mosaic files
    for fitsfile in mosaic_files:
        try:
            beam = read_file(fitsfile).beam
            print(f"Processing mosaic: {os.path.basename(fitsfile)}", flush=True)
            print(f"Beam size: {beam.major.to(u.arcsec):.2f} x {beam.minor.to(u.arcsec):.2f}")
        except Exception as e:
            # warnings.warn(f"Could not get beam from {fitsfile}: {str(e)}")
            continue
            
        for reg in regions:
            flux, peak, rms = get_flux_in_region(fitsfile, reg)
            try:
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
            except Exception as e:
                warnings.warn(f"Error processing region in {fitsfile}: {str(e)}")
    
    mous_pattern = os.path.join(basepath, 'data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001*/calibrated/working/*image.tt0.pbcor.fits')  # Adjust pattern based on MOUS file location
    mous_files = glob.glob(mous_pattern)
    assert len(mous_files) > 40

    print(f"Starting MOUS loop with {len(mous_files)}", flush=True)
    # Process MOUS files
    for fitsfile in mous_files:
        print(f"Processing MOUS: {os.path.basename(fitsfile)}")
        try:
            beam = read_file(fitsfile).beam
            print(f"Beam size: {beam.major.to(u.arcsec):.2f} x {beam.minor.to(u.arcsec):.2f}")
        except Exception as e:
            warnings.warn(f"Could not get beam from {fitsfile}: {str(e)}")
            continue
            
        for reg in regions:
            flux, peak, rms = get_flux_in_region(fitsfile, reg)
            try:
                results.append({
                    'filename': os.path.basename(fitsfile),
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
    output_file = os.path.join(basepath, 'tables', 'point_source_fluxes.ecsv')
    results_table.write(output_file, format='ascii.ecsv', overwrite=True)
    print(f"Results saved to {output_file}")

if __name__ == '__main__':
    main()
