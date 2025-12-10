"""
Reproject final products to RA/Dec.

Target header is header_12m_radec.hdr

target folder is {ACES_rootdir}/products_for_ALMA/group.uid___A001_X1590_X30a9.lp_slongmore/.  All non-PV files in this directory should be reprojected.

renaming: reprojected files should go in {ACES_rootdir}//products_for_ALMA/group.uid___A001_X1590_X30a9.lp_slongmore/icrs/, and .icrs should be added to each file's name after `cmz_mosaic`
"""

import os
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
import numpy as np
import warnings
from tqdm import tqdm

# Configuration
ACES_rootdir = Path("/orange/adamginsburg/ACES")
header_file = ACES_rootdir / "reduction_ACES/aces/imaging/data/header_12m_radec.hdr"
source_dir = ACES_rootdir / "products_for_ALMA/group.uid___A001_X1590_X30a9.lp_slongmore"
output_dir = source_dir / "icrs"


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


def reproject_file(input_path, target_header, output_path):
    """
    Reproject a FITS file to the target coordinate system.
    
    Parameters:
        input_path: Path to input FITS file
        target_header: Target header for reprojection
        output_path: Path for output FITS file
    """
    print(f"\nReprojecting: {input_path.name}")
    
    with fits.open(input_path) as hdul:
        hdu = hdul[0]

        # Reproject the image
        reprojected, footprint = reproject_interp(hdu, target_header, return_footprint=True)
        
        new_header = hdu.header.copy()
        new_header.update(target_header)
        
        new_hdu = fits.PrimaryHDU(data=reprojected, header=new_header)
    
        # Write the output file
        print(f"    Writing to: {output_path.name}")
        new_hdu.writeto(output_path, overwrite=True)
        print(f"    ✓ Complete")


def main():
    """Main reprojection workflow."""
    print("ACES Data Reprojection to ICRS (RA/Dec)")
    
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
        
        reproject_file(input_path, target_header, output_path)
    
    # Summary
    print(f"Reprojection complete!  Total: {len(files_to_reproject)}")

if __name__ == "__main__":
    main()