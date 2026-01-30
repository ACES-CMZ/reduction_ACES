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
from radio_beam import Beam
import radio_beam

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

    # VLA radio continuum - multiple files covering different regions
    'VLA 6cm': [
        '/orange/adamginsburg/cmz/xinglu/20kms_VLA_continuum_IM.fits',
        '/orange/adamginsburg/cmz/xinglu/c20kms_cont_merged.fits',
        '/orange/adamginsburg/cmz/xinglu/c20kms_cont_pbcor_merged.fits',
        '/orange/adamginsburg/cmz/xinglu/c50kms_cont.image.fits',
        '/orange/adamginsburg/cmz/xinglu/c50kms_cont.image.pbcor.fits',
        '/orange/adamginsburg/cmz/xinglu/DustRidge_CONT_tclean_nterm2.image.tt0.fits',
        '/orange/adamginsburg/cmz/xinglu/SgrA_CONT_tclean_nterm2.image.tt0.fits',
        '/orange/adamginsburg/cmz/xinglu/sgrb1off_cont_merged.fits',
        '/orange/adamginsburg/cmz/xinglu/sgrb1off_cont_pbcor_merged.fits',
        '/orange/adamginsburg/cmz/xinglu/SgrB1off_VLA_continuum_IM.fits',
        '/orange/adamginsburg/cmz/xinglu/sgrc_cont_merged.fits',
        '/orange/adamginsburg/cmz/xinglu/sgrc_cont_pbcor_merged.fits',
        '/orange/adamginsburg/cmz/xinglu/SgrC_CONT_tclean_nterm2.pbcor.image.tt0.fits',
    ],
    'MEERKAT 20cm': '/orange/adamginsburg/cmz/meerkat/MeerKAT_Galactic_Centre_1284MHz-StokesI.fits'
}

# Filter to only existing files
existing_files = {}
for label, filepath in data_files.items():
    if isinstance(filepath, list):
        # For lists (e.g., VLA 6cm), keep all existing files
        existing = [f for f in filepath if Path(f).exists()]
        if existing:
            existing_files[label] = existing
        else:
            print(f"Warning: No files found for {label}")
    else:
        if Path(filepath).exists():
            existing_files[label] = filepath
        else:
            print(f"Warning: File not found: {filepath}")

print(f"\nFound {len(existing_files)} existing data files")

# Define wavelengths and beams for SED extraction
wavelengths = {
    'Spitzer 3.6μm': 3.6*u.um,
    'Spitzer 4.5μm': 4.5*u.um,
    'Spitzer 5.8μm': 5.8*u.um,
    'Spitzer 8.0μm': 8.0*u.um,
    '24um': 24*u.um,
    '70um': 70*u.um,
    'Herschel 70um': 70*u.um,
    '160um': 160*u.um,
    'Herschel 160um': 160*u.um,
    '250um': 250*u.um,
    'Herschel 250um': 250*u.um,
    '1mm': (230*u.GHz).to(u.um, u.spectral()),
    'CMZoom 1mm': (230*u.GHz).to(u.um, u.spectral()),
    '3mm': (96*u.GHz).to(u.um, u.spectral()),
    'ACES 3mm': (96*u.GHz).to(u.um, u.spectral()),
    '6cm': 6*u.cm,
    'VLA 6cm': 6*u.cm,
    '20cm': 20*u.cm,
    'VLA 20cm': 20*u.cm,
    'MEERKAT 20cm': 20*u.cm
}

# Beam sizes from Traficante+ 2011 and instrument specs
beam_sizes = {
    'Spitzer 3.6μm': Beam(2*u.arcsec),
    'Spitzer 4.5μm': Beam(2*u.arcsec),
    'Spitzer 5.8μm': Beam(2*u.arcsec),
    'Spitzer 8.0μm': Beam(2*u.arcsec),
    '24um': Beam(6*u.arcsec),
    '70um': Beam(10.7*u.arcsec, 9.7*u.arcsec),
    'Herschel 70um': Beam(10.7*u.arcsec, 9.7*u.arcsec),
    '160um': Beam(13.9*u.arcsec, 13.2*u.arcsec),
    'Herschel 160um': Beam(13.9*u.arcsec, 13.2*u.arcsec),
    '250um': Beam(23.9*u.arcsec, 22.8*u.arcsec),
    'Herschel 250um': Beam(23.9*u.arcsec, 22.8*u.arcsec),
    '1mm': Beam(3*u.arcsec),  # CMZoom typical
    'CMZoom 1mm': Beam(3*u.arcsec),
    '3mm': Beam(1.5*u.arcsec),  # ACES typical
    'ACES 3mm': Beam(1.5*u.arcsec),
    '6cm': None,  # Will be read from FITS header
    'VLA 6cm': None,
    '20cm': None,  # Will be read from FITS header
    'VLA 20cm': None,
    'MEERKAT 20cm': None
}

# Cache for file WCS and bounds to avoid repeated file opening
_file_cache = {}

def select_best_file(coord, filepaths):
    """
    Select the best file from a list based on overlap and beam size.

    Parameters
    ----------
    coord : SkyCoord
        Source coordinate
    filepaths : list
        List of file paths to check

    Returns
    -------
    best_file : str or None
        Path to the best file, or None if no overlap
    """
    candidates = []

    for filepath in filepaths:
        # Check cache first
        if filepath not in _file_cache:
            try:
                with fits.open(filepath) as hdul:
                    # Find image HDU
                    image_hdu = None
                    for h in hdul:
                        if h.data is not None and len(h.data.shape) >= 2:
                            image_hdu = h
                            break

                    if image_hdu is None:
                        _file_cache[filepath] = None
                        continue

                    # Get WCS
                    try:
                        wcs = WCS(image_hdu.header).celestial
                    except (AttributeError, ValueError):
                        wcs = WCS(image_hdu.header)
                        if wcs.naxis > 2:
                            wcs = wcs.sub(['longitude', 'latitude'])

                    # Get data shape
                    data_shape = image_hdu.data.shape
                    while len(data_shape) > 2:
                        data_shape = data_shape[1:]

                    # Get beam size
                    try:
                        beam = radio_beam.Beam.from_fits_header(image_hdu.header)
                        beam_area = beam.sr.value
                    except (radio_beam.NoBeamException, KeyError):
                        # Should not happen for VLA data, but handle gracefully
                        beam_area = 1e10  # Large value = low priority

                    _file_cache[filepath] = (wcs, data_shape, beam_area)
            except Exception:
                # Skip files that can't be opened
                _file_cache[filepath] = None
                continue

        # Use cached info
        cache_info = _file_cache.get(filepath)
        if cache_info is None:
            continue

        wcs, data_shape, beam_area = cache_info

        # Check if coordinate is within image bounds
        try:
            xx, yy = wcs.world_to_pixel(coord)
            if 0 <= xx < data_shape[1] and 0 <= yy < data_shape[0]:
                candidates.append((filepath, beam_area))
        except Exception:
            continue

    if not candidates:
        return None

    # Select file with smallest beam (best resolution)
    candidates.sort(key=lambda x: x[1])
    return candidates[0][0]


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
    except (ValueError, IndexError, AttributeError) as e:
        # Expected: coordinate outside image bounds, invalid WCS, attribute errors
        print(f"    Error creating cutout mask: {type(e).__name__}: {e}")
        return None, None

    # Reproject to Galactic frame with optimal WCS
    try:
        csys, sz = reproject.mosaicking.find_optimal_celestial_wcs([(co, wcs[slcs])], frame='galactic')
        newdata, _ = reproject.reproject_interp(input_data=(co, wcs[slcs]),
                                                output_projection=csys,
                                                shape_out=sz)
    except (ValueError, TypeError, AttributeError) as e:
        # Expected: WCS problems, incompatible projections, shape issues
        print(f"    Error reprojecting to Galactic: {type(e).__name__}: {e}")
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
    except (ValueError, IndexError, AttributeError) as e:
        # Expected: slicing issues, WCS problems
        print(f"    Error extracting final region: {type(e).__name__}: {e}")
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


def extract_sed(coord, data_files, source, detection_sigma=5.0):
    """
    Extract SED by measuring flux at source position in each image.

    Parameters
    ----------
    coord : SkyCoord
        Source coordinate
    data_files : dict
        Dictionary of wavelength labels to file paths
    source : Table row
        Source catalog entry
    detection_sigma : float
        S/N threshold for detections

    Returns
    -------
    sed_table : dict
        Dictionary with wavelengths, fluxes, errors, and detection flags
    """
    sed_data = {'wavelength': [], 'flux': [], 'error': [], 'is_detection': [], 'label': []}

    # Always add ACES 3mm point from catalog (this is a detection by definition)
    if 'flux' in source.colnames and source['flux'] > 0:
        sed_data['wavelength'].append((96*u.GHz).to(u.um, u.spectral()))
        sed_data['flux'].append(source['flux'] * u.Jy)
        sed_data['error'].append(source.get('flux_err', source['flux']*0.1) * u.Jy)
        sed_data['is_detection'].append(True)
        sed_data['label'].append('3mm (ACES)')

    # Add CMZoom 1mm if available
    if 'cmzoom_flux' in source.colnames and not np.isnan(source['cmzoom_flux']) and source['cmzoom_flux'] > 0:
        cmz_err = source.get('cmzoom_flux_err', source['cmzoom_flux']*0.2)
        is_det = (source['cmzoom_flux'] / cmz_err) > detection_sigma if not np.isnan(cmz_err) else False
        sed_data['wavelength'].append((230*u.GHz).to(u.um, u.spectral()))
        sed_data['flux'].append(source['cmzoom_flux'] * u.Jy)
        sed_data['error'].append(cmz_err * u.Jy)
        sed_data['is_detection'].append(is_det)
        sed_data['label'].append('1mm (CMZoom)')

    # Extract from other data files
    for label, filepath in data_files.items():
        # Skip if we already added from catalog
        if label in ['3mm', '1mm', 'ACES 3mm', 'CMZoom 1mm']:
            continue

        # Handle lists of files (e.g., VLA 6cm)
        if isinstance(filepath, list):
            best_file = select_best_file(coord, filepath)
            if best_file is None:
                print(f"    SED: {label} - no overlapping files")
                continue
            filepath = best_file

        try:
            with fits.open(filepath) as hdul:
                # Find image HDU
                image_hdu = None
                for h in hdul:
                    if h.data is not None and len(h.data.shape) >= 2:
                        image_hdu = h
                        break

                if image_hdu is None:
                    print(f"    SED: No valid HDU in {label}")
                    continue

                # Get data and WCS
                data = image_hdu.data.copy()
                header = image_hdu.header

                # Handle multi-dimensional data
                while len(data.shape) > 2:
                    data = data[0]

                # Get WCS
                try:
                    wcs = WCS(header).celestial
                except (AttributeError, ValueError) as e:
                    # If celestial WCS extraction fails, try subsetting
                    print(f"    SED: {label} - WCS.celestial failed ({type(e).__name__}), trying sub")
                    wcs = WCS(header)
                    if wcs.naxis > 2:
                        wcs = wcs.sub(['longitude', 'latitude'])

                # Get pixel coordinates
                xx, yy = wcs.world_to_pixel(coord)
                xx, yy = int(np.round(xx)), int(np.round(yy))

                # Check bounds
                if not (0 <= xx < data.shape[1] and 0 <= yy < data.shape[0]):
                    print(f"    SED: {label} - coordinate outside image bounds")
                    continue

                val = data[yy, xx]

                if np.isnan(val) or not np.isfinite(val):
                    print(f"    SED: {label} - NaN or non-finite value at position")
                    continue

                # Get beam for flux conversion
                beam = beam_sizes.get(label)
                if beam is None:
                    if '6cm' in label or '20cm' in label:
                        try:
                            beam = radio_beam.Beam.from_fits_header(header)
                        except (radio_beam.NoBeamException, KeyError) as e:
                            # VLA/MeerKAT data should always have beams - this indicates a real problem
                            print(f"    SED: {label} - ERROR: missing beam in header ({type(e).__name__}: {e})")
                            # this outcome isn't OK.
                            raise e
                    else:
                        print(f"    SED: {label} - no beam size defined, skipping")
                        continue

                # Convert to Jy based on wavelength and units
                bunit = str(header.get('BUNIT', '')).upper()

                if wavelengths[label] < 1*u.mm:  # IR: typically in MJy/sr
                    if 'MJY/SR' in bunit or 'MJY.SR' in bunit:
                        flux = (val * u.MJy/u.sr * beam.sr).to(u.Jy)
                    elif 'JY' in bunit:
                        flux = val * u.Jy
                    else:
                        # Assume MJy/sr for IR
                        flux = (val * u.MJy/u.sr * beam.sr).to(u.Jy)
                else:  # Radio: typically in Jy/beam or Jy
                    if 'JY/BEAM' in bunit:
                        flux = val * u.Jy
                    elif 'JY' in bunit:
                        flux = val * u.Jy
                    else:
                        # Assume Jy/beam for radio
                        flux = val * u.Jy

                # For upper limits, use absolute value
                flux_val = abs(flux.value) if flux.value < 0 else flux.value
                flux = flux_val * u.Jy

                # Check if this is a significant detection at the catalog position
                # Calculate SNR using local background subtraction
                is_detection = False
                peak_snr = 0
                try:
                    # Extract local region for background estimation
                    # Use annulus: inner radius excludes source (3 pixels), outer radius for background
                    half_size = 10  # pixels for background region
                    x_min = max(0, xx - half_size)
                    x_max = min(data.shape[1], xx + half_size + 1)
                    y_min = max(0, yy - half_size)
                    y_max = min(data.shape[0], yy + half_size + 1)

                    local_region = data[y_min:y_max, x_min:x_max]

                    if np.any(np.isfinite(local_region)):
                        # Create mask excluding central 3x3 pixels (inner aperture containing source)
                        center_y = yy - y_min
                        center_x = xx - x_min
                        bg_mask = np.ones_like(local_region, dtype=bool)

                        # Exclude inner 3x3 region
                        for dy in range(-1, 2):
                            for dx in range(-1, 2):
                                cy, cx = center_y + dy, center_x + dx
                                if 0 <= cy < bg_mask.shape[0] and 0 <= cx < bg_mask.shape[1]:
                                    bg_mask[cy, cx] = False

                        # Get background pixels
                        bg_pixels = local_region[bg_mask & np.isfinite(local_region)]

                        if len(bg_pixels) > 10:
                            # Calculate background statistics using median and MAD
                            bg_median = np.median(bg_pixels)
                            mad = np.median(np.abs(bg_pixels - bg_median))
                            bg_rms = 1.4826 * mad  # Convert MAD to standard deviation

                            # Calculate SNR at the catalog position
                            if bg_rms > 0:
                                snr = (val - bg_median) / bg_rms
                                peak_snr = snr

                                # Detection if SNR > 3
                                if snr > 3.0:
                                    is_detection = True
                except (IndexError, ValueError, TypeError, AttributeError) as e:
                    # If detection check fails, assume upper limit
                    print(f"    SED: {label} - detection check failed ({type(e).__name__})")
                    is_detection = False

                # Estimate error based on detection status
                if is_detection:
                    # For detections, use 20% calibration uncertainty
                    error = max(flux * 0.2, 0.001*u.Jy)
                else:
                    # For upper limits, use 30% uncertainty
                    error = max(flux * 0.3, 0.001*u.Jy)

                sed_data['wavelength'].append(wavelengths[label])
                sed_data['flux'].append(flux)
                sed_data['error'].append(error)
                sed_data['is_detection'].append(is_detection)
                sed_data['label'].append(label)

                det_type = f"DETECTION (SNR={peak_snr:.1f})" if is_detection else "upper limit"
                print(f"    SED: {label} - flux={flux.value:.4g} Jy, det={is_detection} ({det_type})")

        except (OSError, IOError, ValueError, KeyError) as e:
            # Expected errors: file issues, WCS problems, unit conversions
            print(f"    SED: Error processing {label}: {type(e).__name__}: {e}")
            continue

    return sed_data


def modified_blackbody_simple(wavelength, temperature, normalization):
    """
    Simple modified blackbody with fixed beta=1.5.

    Parameters
    ----------
    wavelength : Quantity
        Wavelengths to evaluate
    temperature : Quantity
        Dust temperature
    normalization : Quantity
        Flux at reference wavelength to normalize to

    Returns
    -------
    flux : Quantity
        Model flux at each wavelength
    """
    from astropy.modeling.models import BlackBody

    nu = wavelength.to(u.Hz, u.spectral())
    nu_ref = (3*u.mm).to(u.Hz, u.spectral())

    beta = 1.5

    bb = BlackBody(temperature)
    bb_ratio = bb(nu) / bb(nu_ref)

    # Include dust opacity scaling
    opacity_ratio = (nu / nu_ref) ** beta

    flux = normalization * bb_ratio * opacity_ratio

    return flux


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
        # Handle lists of files (e.g., VLA 6cm)
        if isinstance(filepath, list):
            best_file = select_best_file(coord, filepath)
            if best_file is None:
                print(f"  Skipping {label}: no overlapping files")
                continue
            filepath = best_file
            print(f"  Using {Path(filepath).name} for {label}")

        try:
            with fits.open(filepath) as hdul:
                cutout_data, cutout_wcs = galactic_cutout(coord, hdul, size=size)

                # Only include if cutout has valid data
                if cutout_data is not None and not np.all(np.isnan(cutout_data)):
                    cutouts[label] = (cutout_data, cutout_wcs)
                else:
                    print(f"  Skipping {label}: no valid data")
        except (OSError, IOError, ValueError) as e:
            # Expected errors: file issues, WCS/reprojection problems
            print(f"  Skipping {label}: {type(e).__name__}: {e}")

    if len(cutouts) == 0:
        print(f"  No valid cutouts for source {source['index']}, skipping")
        return

    # Extract SED
    sed_data = extract_sed(coord, data_files, source)

    # Create figure with gridspec (reserve bottom-right for SED)
    n_panels = len(cutouts)
    n_cols = 4
    n_rows = int(np.ceil((n_panels + 1) / n_cols))  # +1 for SED panel

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

    # Add SED panel in bottom-right position
    if len(sed_data['wavelength']) > 0:
        # Determine position for SED (bottom-right corner)
        sed_row = n_rows - 1
        sed_col = n_cols - 1

        ax_sed = fig.add_subplot(gs[sed_row, sed_col])

        # Sort by wavelength for plotting
        sorted_idx = np.argsort([w.value for w in sed_data['wavelength']])
        wl_sorted = [sed_data['wavelength'][i] for i in sorted_idx]
        flux_sorted = [sed_data['flux'][i] for i in sorted_idx]
        error_sorted = [sed_data['error'][i] for i in sorted_idx]
        is_det_sorted = [sed_data['is_detection'][i] for i in sorted_idx]

        # Convert to arrays for plotting
        wl_vals = u.Quantity(wl_sorted)
        flux_vals = u.Quantity(flux_sorted)
        error_vals = u.Quantity(error_sorted)

        # Plot detections and upper limits
        for i, (wl, flux, err, is_det) in enumerate(zip(wl_vals, flux_vals, error_vals, is_det_sorted)):
            if is_det:
                # Detection: plot as square
                ax_sed.errorbar(wl.to(u.um).value, flux.to(u.Jy).value,
                              yerr=err.to(u.Jy).value,
                              marker='s', markersize=6, color='blue',
                              markeredgecolor='k', capsize=3, capthick=1.5)
            else:
                # Upper limit: plot as downward arrow
                ax_sed.errorbar(wl.to(u.um).value, flux.to(u.Jy).value,
                              yerr=err.to(u.Jy).value, uplims=True,
                              marker='v', markersize=6, color='gray',
                              markerfacecolor='none', markeredgecolor='k',
                              capsize=3, capthick=1.5)

        # Add 50K modified blackbody through ACES 3mm point
        if 'flux' in source.colnames and source['flux'] > 0:
            model_wl = np.geomspace(1*u.um, 1*u.m, 500)
            model_flux = modified_blackbody_simple(
                model_wl,
                temperature=50*u.K,
                normalization=source['flux']*u.Jy
            )
            ax_sed.plot(model_wl.to(u.um).value, model_flux.to(u.Jy).value,
                       'r--', linewidth=1.5, alpha=0.7, label='50K MBB', zorder=-5)

        ax_sed.set_xscale('log')
        ax_sed.set_yscale('log')
        ax_sed.set_xlabel('Wavelength [μm]', fontsize=11)
        ax_sed.set_ylabel('Flux Density [Jy]', fontsize=11)
        ax_sed.set_title('SED', fontsize=12, fontweight='bold')
        ax_sed.grid(True, alpha=0.3)
        ax_sed.legend(loc='best', fontsize=9)

        # Set reasonable axis limits
        if len(wl_vals) > 0:
            ax_sed.set_xlim(wl_vals.min().to(u.um).value * 0.5,
                           wl_vals.max().to(u.um).value * 2)
            flux_min = flux_vals[flux_vals > 0].min().to(u.Jy).value if np.any(flux_vals > 0) else 1e-3
            flux_max = flux_vals.max().to(u.Jy).value
            ax_sed.set_ylim(flux_min * 0.3, flux_max * 3)

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
        except (KeyboardInterrupt, SystemExit):
            # Don't catch user interrupts
            raise
        except Exception as e:
            # Catch any unexpected errors to continue processing
            print(f"ERROR: Unexpected error processing source {i}: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            continue

    print(f"\nDone! Created cutouts in {output_dir}")


if __name__ == '__main__':
    main()
