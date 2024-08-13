from astropy.io import fits
import numpy as np
import os
import sys
import glob
import pvextractor
from pvextractor import Path, extract_pv_slice
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
import spectral_cube
from spectral_cube import SpectralCube

save_path = '/orange/adamginsburg/ACES/broadline_sources/EVFs/images/'

# Latitude
B_MIN = -0.27*u.deg
B_MAX = 0.22*u.deg
LIST_B = make_position_list(B_MIN, B_MAX, 0.5*u.arcmin, u.arcmin)

# Longitude
L_MIN = -0.59*u.deg
L_MAX = 0.88*u.deg
LIST_L = make_position_list(L_MIN, L_MAX, 0.5*u.arcmin, u.arcmin)

def make_position_list(amin, amax, step, unit):
    return (np.arange(amin.to(unit).value, amax.to(unit).value+1, step.value)*unit).to(u.deg)

def make_pv_b(cube, b, mol):
    reg = regions.RectangleSkyRegion(center=SkyCoord((L_MIN+L_MAX)/2., b, frame='galactic'), width=1.5*u.deg, height=1*u.arcmin)
    subcube = cube.subcube_from_regions([reg])
    pv_mean = subcube.mean(axis=1)
    pv_mean.save(f'{save_path}/{mol}_pv_b{round(b.value, 3)}_mean.fits')
    pv_max = subcube.max(axis=1)
    pv_max.save(f'{save_path}/{mol}_pv_b{round(b.value, 3)}_max.fits')

def make_pv_l(cube, l, mol):
    reg = regions.RectangleSkyRegion(center=SkyCoord(l, (B_MIN+B_MAX)/2., frame='galactic'), width=1*u.arcmin, height=0.5*u.deg)
    subcube = cube.subcube_from_regions([reg])
    pv_mean = subcube.mean(axis=2)
    pv_mean.save(f'{save_path}/{mol}_pv_l{round(l.value, 3)}_mean.fits')
    pv_max = subcube.max(axis=2)
    pv_max.save(f'{save_path}/{mol}_pv_l{round(l.value, 3)}_max.fits')

def make_pv_mol(cube_fn):
    cube = SpectralCube.read(cube_fn)
    mol = cube_fn.split('/')[-1].split('_')[0]

    for b in list_b:
        make_pv_b(cube, b, mol)

    for l in list_l:
        make_pv_l(cube, l, mol)

def main():
    basepath = '/orange/adamginsburg/ACES/mosaics/cubes/'
    cube_fn_CS = f'{basepath}/CS_CubeMosaic.fits'
    
    make_pv_mol(cube_fn_CS)


if __name__ == "__main__":
    main()#sys.argv[1])