from astropy.io import fits
from astropy.wcs import WCS
import itertools
import os
import warnings
import numpy as np

import aces
basepath = aces.conf.basepath

molecules = ['CS21', 'HC3N', 'HCOP', 'HNCO_7m12mTP', 'SiO21', 'SO21', 'SO32']

# Ratios with a specific numerator/denominator direction that the pairwise
# combinations above do not produce (combinations use molecules-list order, so
# they yield CS21_over_HC3N etc.; these give the inverse direction requested).
extra_pairs = [('HNCO_7m12mTP', 'CS21'), ('HC3N', 'CS21')]

moments_dir = f'{basepath}/mosaics/cubes/moments'
cubes_dir = f'{basepath}/mosaics/cubes'


def compute_downsampled_moments(mol):
    """Make max + mom0 moment maps from a molecule's spatially-downsampled cube
    (``*_CubeMosaic_downsampled9.fits``), so ratios can be built at the smoothed
    resolution.  No-op if the downsampled cube is absent or the moments exist."""
    cube_fn = f'{cubes_dir}/{mol}_CubeMosaic_downsampled9.fits'
    if not os.path.exists(cube_fn):
        return
    hdr = fits.getheader(cube_fn)
    cd3 = abs(hdr['CDELT3'])
    if 'm/s' in str(hdr.get('CUNIT3', '')).lower():
        cd3 /= 1000.0  # -> km/s
    h2 = WCS(hdr).celestial.to_header()
    cube = None
    for suf in ('max', 'mom0'):
        out = f'{moments_dir}/{mol}_CubeMosaic_downsampled9_{suf}.fits'
        if os.path.exists(out):
            continue
        if cube is None:
            cube = fits.getdata(cube_fn)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            m = np.nanmax(cube, axis=0) if suf == 'max' else np.nansum(cube, axis=0) * cd3
        fits.PrimaryHDU(data=m, header=h2).writeto(out, overwrite=True)
        print(f"wrote downsampled moment {mol}_downsampled9_{suf}")


def make_ratio(mol1, mol2, suffix):
    f1 = f'{moments_dir}/{mol1}_CubeMosaic_{suffix}.fits'
    f2 = f'{moments_dir}/{mol2}_CubeMosaic_{suffix}.fits'
    if not (os.path.exists(f1) and os.path.exists(f2)):
        print(f"skip {mol1}_over_{mol2}_{suffix}: missing input")
        return
    fh1 = fits.open(f1)
    fh2 = fits.open(f2)
    if np.all(fh2[0].data.shape == fh1[0].data.shape):
        fh1[0].data = fh1[0].data / fh2[0].data
        out = f'{moments_dir}/ratios/{mol1}_over_{mol2}_{suffix}.fits'
        fh1.writeto(out, overwrite=True)
        print(f"wrote {mol1}_over_{mol2}_{suffix}")


if __name__ == '__main__':
    os.makedirs(f'{moments_dir}/ratios', exist_ok=True)
    # build the smoothed-resolution moment maps first
    for mol in molecules:
        compute_downsampled_moments(mol)
    # ratios at full resolution and at the downsampled (smoothed) resolution
    for suffix in ('max', 'mom0', 'downsampled9_max', 'downsampled9_mom0'):
        for mol1, mol2 in itertools.combinations(molecules, 2):
            make_ratio(mol1, mol2, suffix)
        for mol1, mol2 in extra_pairs:
            make_ratio(mol1, mol2, suffix)
