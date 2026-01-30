"""
The "final" catalog of sources in ACES is here:
/orange/adamginsburg/ACES/upload/compact_cont_source_catalog/official_cont_catalog_files/ACES_catalog_v0_README.txt
/orange/adamginsburg/ACES/upload/compact_cont_source_catalog/official_cont_catalog_files/aces_compact_catalog_v0.fits

Following the example of /orange/adamginsburg/ACES/broadline_sources/G0.025-0.073/MUBLO_MultiwavelengthCutouts.ipynb, create cutout visualizations for each source in the catalog.

some of the file paths may be different.  For example, the Spitzer images can be found at /orange/adamginsburg/galactic_plane_surveys/glimpse/, and instead of the 850um data, we can use the CMZoom data: /orange/adamginsburg/cmz/cmzoom/continuum_mosaic/CMZoom_continuum_pbcor.fits

For each of the ancillary data sets, the panel should only be created if the selected source is covered; if the cutout is blank or all NaN, it should not be included.

The source should be marked by an ellipse on each panel, with the major and minor axes and position angle taken from the catalog.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy import units as u
from reproject import reproject_interp
import regions
import warnings
from pathlib import Path

# Suppress warnings
warnings.filterwarnings('ignore')

# Define output directory
output_dir = Path('/orange/adamginsburg/ACES/figures/cutouts')
output_dir.mkdir(parents=True, exist_ok=True)

# Load catalog
catalog_path = '/orange/adamginsburg/ACES/tables/aces_compact_catalog_v0_withExtras.fits'
catalog = Table.read(catalog_path)
print(f"Loaded {len(catalog)} sources from catalog")

# Define data files for different wavelengths
# Using exact paths from MUBLO notebook
data_files = {
    # Spitzer IRAC bands (3.6, 4.5, 5.8, 8.0 μm)
    'Spitzer 3.6μm': '/orange/adamginsburg/galactic_plane_surveys/glimpse/GLM_00000+0000_mosaic_I1.fits',
    'Spitzer 4.5μm': '/orange/adamginsburg/galactic_plane_surveys/glimpse/GLM_00000+0000_mosaic_I2.fits',
    'Spitzer 5.8μm': '/orange/adamginsburg/galactic_plane_surveys/glimpse/GLM_00000+0000_mosaic_I3.fits',
    'Spitzer 8.0μm': '/orange/adamginsburg/galactic_plane_surveys/glimpse/GLM_00000+0000_mosaic_I4.fits',

    # alt path /orange/adamginsburg/spitzer/GLIMPSE/GLM_00000+0000_mosaic_I1.fits
    # alt path /orange/adamginsburg/spitzer/GLIMPSE/GLM_00000+0000_mosaic_I2.fits
    # alt path /orange/adamginsburg/spitzer/GLIMPSE/GLM_00000+0000_mosaic_I3.fits
    # alt path /orange/adamginsburg/spitzer/GLIMPSE/GLM_00000+0000_mosaic_I4.fits
    
    # MIPS 24μm
    '24um': '/orange/adamginsburg/cmz/mipsgal_24micron_data/gc_mosaic_MIPSGAL_gal.fits',
    
    # Herschel PACS and SPIRE bands
    'Herschel 70um': '/orange/adamginsburg/cmz/herschel/destripe_l000_blue_wgls_rcal.fits',
    'Herschel 160um': '/orange/adamginsburg/cmz/herschel/destripe_l000_red_wgls_rcal.fits',
    'Herschel 250um': '/orange/adamginsburg/cmz/herschel/destripe_l000_PSW_wgls_rcal.fits',
    
    # CMZoom 1mm continuum
    'CMZoom 1mm': '/orange/adamginsburg/cmz/cmzoom/continuum_mosaic/CMZoom_continuum_pbcor.fits',
    
    # ACES 3mm continuum
    'ACES 3mm': '/orange/adamginsburg/ACES/mosaics/continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits',
    
    # VLA radio continuum
    # for VLA 6cm, see also: /orange/adamginsburg/cmz/xinglu/*fits
    'VLA 6cm': '/orange/adamginsburg/cmz/xinglu/SgrA_CONT_tclean_nterm2.image.tt0.fits',
    'VLA 20cm': '/orange/adamginsburg/cmz/meerkat/MeerKAT_Galactic_Centre_1284MHz-StokesI.fits'
}

# Filter to only existing files
existing_files = {}
for label, filepath in data_files.items():
    if Path(filepath).exists():
        existing_files[label] = filepath
    else:
        print(f"Warning: File not found: {filepath}")

print(f"\nFound {len(existing_files)} existing data files")

def galactic_cutout(coord, hdu, size=50*u.arcsec):
    """
    Extract a cutout from a FITS file and reproject to Galactic coordinates.
    Based on the MUBLO notebook implementation.
    
    Parameters
    ----------
    coord : SkyCoord
        Coordinate of the source
    hdu : astropy.io.fits.ImageHDU or HDUList
        FITS HDU containing the image and WCS
    size : Quantity
        Size of the cutout (width/height of square region)
        
    Returns
    -------
    cutout_data : ndarray
        Reprojected cutout data
    cutout_wcs : WCS
        WCS of the cutout in Galactic coordinates
    """
    import reproject.mosaicking
    import regions
    
    # Get the image HDU
    if isinstance(hdu, fits.HDUList):
        # Find first valid image HDU
        image_hdu = None
        for h in hdu:
            if h.data is not None and len(h.data.shape) >= 2:
                image_hdu = h
                break
        if image_hdu is None:
            return None, None
    else:
        image_hdu = hdu
    
    # Handle FITS cubes (take first slice)
    data = image_hdu.data
    wcs = WCS(image_hdu.header).celestial
    
    # Squeeze data to 2D
    while len(data.shape) > 2:
        data = data[0]
    
    # Create a circular region for the cutout (circumscribing the target square)
    cutregion = regions.CircleSkyRegion(center=coord, radius=size * np.sqrt(2) / 2)
    
    # Make cutout in original projection
    try:
        msk = cutregion.to_pixel(wcs).to_mask()
        slcs, _ = msk.get_overlap_slices(data.shape)
        if slcs is None:
            return None, None
        co = msk.cutout(data)
        if co is None:
            return None, None
    except Exception as e:
        print(f"    Error creating cutout mask: {e}")
        return None, None
    
    # Reproject to Galactic frame with optimal WCS
    try:
        csys, sz = reproject.mosaicking.find_optimal_celestial_wcs([(co, wcs[slcs])], frame='galactic')
        newdata, _ = reproject.reproject_interp(input_data=(co, wcs[slcs]),
                                                output_projection=csys,
                                                shape_out=sz)
    except Exception as e:
        print(f"    Error reprojecting to Galactic: {e}")
        return None, None
    
    # Create the target rectangular region in the reprojected data
    region = regions.RectangleSkyRegion(center=coord, width=size, height=size)
    try:
        preg = region.to_pixel(csys)
        msk = preg.to_mask()
        slcs2, _ = msk.get_overlap_slices(newdata.shape)
        if slcs2 is None:
            return None, None
        final_data = newdata[slcs2]
        final_wcs = csys[slcs2]
    except Exception as e:
        print(f"    Error extracting final region: {e}")
        return None, None
    
    return final_data, final_wcs


def create_custom_colormap():
    """
    Create a custom colormap combining gray_r for low values and hot for high values.
    """
    # Get base colormaps
    gray_r = plt.cm.gray_r
    hot = plt.cm.hot
    
    # Sample them
    n_samples = 128
    colors_gray = gray_r(np.linspace(0, 1, n_samples))
    colors_hot = hot(np.linspace(0, 1, n_samples))
    
    # Combine: gray_r for lower half, hot for upper half
    colors = np.vstack([colors_gray, colors_hot])
    
    # Create new colormap
    custom_cmap = LinearSegmentedColormap.from_list('gray_hot', colors)
    
    return custom_cmap


def plot_multiwavelength_cutout(source_idx, catalog, data_files, output_dir, 
                                 size=50*u.arcsec, save=True):
    """
    Create a multiwavelength cutout visualization for a single source.
    
    Parameters
    ----------
    source_idx : int
        Index of the source in the catalog
    catalog : Table
        ACES catalog
    data_files : dict
        Dictionary mapping labels to file paths
    output_dir : Path
        Output directory for saved figures
    size : Quantity
        Size of cutouts
    save : bool
        Whether to save the figure
    """
    source = catalog[source_idx]
    coord = SkyCoord(source['GLON_peak'], source['GLAT_peak'], 
                     unit=u.deg, frame='galactic')
    
    print(f"\nProcessing source {source['index']} at ({source['GLON_peak']:.4f}, {source['GLAT_peak']:.4f})")
    
    # Prepare cutouts
    cutouts = {}
    for label, filepath in data_files.items():
        try:
            with fits.open(filepath) as hdul:
                cutout_data, cutout_wcs = galactic_cutout(coord, hdul, size=size)
                
                # Only include if cutout has valid data
                if cutout_data is not None and not np.all(np.isnan(cutout_data)):
                    cutouts[label] = (cutout_data, cutout_wcs)
                else:
                    print(f"  Skipping {label}: no valid data")
        except Exception as e:
            print(f"  Skipping {label}: {e}")
    
    if len(cutouts) == 0:
        print(f"  No valid cutouts for source {source['index']}, skipping")
        return
    
    # Create figure with gridspec
    n_panels = len(cutouts)
    n_cols = 4
    n_rows = int(np.ceil(n_panels / n_cols))
    
    fig = plt.figure(figsize=(16, 4 * n_rows))
    gs = gridspec.GridSpec(n_rows, n_cols, figure=fig, hspace=0.3, wspace=0.3)
    
    # Custom colormap
    cmap = create_custom_colormap()
    
    # Plot each cutout
    for i, (label, (data, wcs)) in enumerate(cutouts.items()):
        row = i // n_cols
        col = i % n_cols
        
        ax = fig.add_subplot(gs[row, col], projection=wcs)
        
        # Compute percentile-based scaling for better visualization
        vmin = np.nanpercentile(data, 1)
        vmax = np.nanpercentile(data, 99.5)
        
        # Plot the image
        im = ax.imshow(data, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
        
        # Mark the source with scatter_coord (like MUBLO notebook)
        ax.scatter_coord(coord, marker='o', facecolor='none', edgecolor='g', s=100)
        
        # Also add ellipse if we have fitted parameters
        if 'fitted_major' in source.colnames and not np.isnan(source['fitted_major']):
            from matplotlib.patches import Ellipse
            
            # Get ellipse parameters
            major = source['fitted_major'] * u.arcsec
            minor = source['fitted_minor'] * u.arcsec
            pa = source['pa'] * u.deg
            
            # Convert to pixel coordinates for the ellipse center
            center_pix = wcs.world_to_pixel(coord)
            
            # Convert sizes to pixels
            pixel_scale = np.abs(wcs.wcs.cdelt[0]) * u.deg
            major_pix = (major / pixel_scale).decompose().value
            minor_pix = (minor / pixel_scale).decompose().value
            
            # Create ellipse (PA measured from North through East in Galactic)
            ellipse = Ellipse(xy=center_pix, width=major_pix, height=minor_pix,
                            angle=pa.value, edgecolor='cyan', facecolor='none',
                            linewidth=1.5, linestyle='--', alpha=0.7)
            ax.add_patch(ellipse)
        
        # Set labels
        ax.set_xlabel('Galactic Longitude')
        ax.set_ylabel('Galactic Latitude')
        ax.set_title(label, fontsize=12, fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Intensity', fontsize=10)
    
    # Add overall title
    title = f"Source {source['index']}: GLON={source['GLON_peak']:.4f}°, GLAT={source['GLAT_peak']:.4f}°"
    if 'flux' in source.colnames:
        title += f"\nFlux (3mm) = {source['flux']:.3f} Jy"
    if 'cmzoom_flux' in source.colnames and not np.isnan(source['cmzoom_flux']):
        title += f", Flux (1mm) = {source['cmzoom_flux']:.3f} Jy"
    
    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
    
    if save:
        output_file = output_dir / f"source_{source['index']:04d}_cutouts.png"
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"  Saved to {output_file}")
        plt.close(fig)
    else:
        plt.show()


def main():
    """
    Main function to create cutouts for all sources.
    """
    print(f"Creating cutouts for {len(catalog)} sources")
    print(f"Output directory: {output_dir}")
    
    # Process each source
    for i in range(len(catalog)):
        try:
            plot_multiwavelength_cutout(i, catalog, existing_files, output_dir,
                                       size=50*u.arcsec, save=True)
        except Exception as e:
            print(f"Error processing source {i}: {e}")
            continue
    
    print(f"\nDone! Created cutouts in {output_dir}")


if __name__ == '__main__':
    main()
