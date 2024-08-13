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

# Path to where to save the PV diagrams
save_path = '/orange/adamginsburg/ACES/broadline_sources/EVFs/images/'

# Latitude extrema
B_MIN = -0.27*u.deg
B_MAX = 0.22*u.deg

# Longitude extrema
L_MIN = -0.59*u.deg
L_MAX = 0.88*u.deg

def make_position_list(amin, amax, step=0.5*u.arcmin, unit=u.arcmin):
    """ 
    Make a list of Sky positions from amin to amax with a given step size and unit.
    Ideally used to make a list of latitudes or longitudes for the centers of region cutouts.

    Parameters
    ----------
    amin : Quantity
        Minimum value of the list.
    amax : Quantity
        Maximum value of the list.
    step : Quantity
        Step size between values.
    unit : astropy.unit
        Unit of the values.
    """

    return (np.arange(amin.to(unit).value, amax.to(unit).value+1, step.value)*unit).to(u.deg)

def make_pv_b(cube, b, mol):
    """ 
    Make a PV diagram along a given latitude.

    Parameters
    ----------
    cube : SpectralCube
        SpectralCube object of the data.
    b : Quantity
        Latitude of the center of the region.
    mol : str
        Molecule name.
    """

    reg = regions.RectangleSkyRegion(center=SkyCoord((L_MIN+L_MAX)/2., b, frame='galactic'), width=1.5*u.deg, height=1*u.arcmin)
    subcube = cube.subcube_from_regions([reg])
    pv_mean = subcube.mean(axis=1)
    pv_mean.save(f'{save_path}/{mol}_pv_b{round(b.value, 3)}_mean.fits')
    pv_max = subcube.max(axis=1)
    pv_max.save(f'{save_path}/{mol}_pv_b{round(b.value, 3)}_max.fits')

def make_pv_l(cube, l, mol):
    """ 
    Make a PV diagram along a given longitude.

    Parameters
    ----------
    cube : SpectralCube
        SpectralCube object of the data.
    l : Quantity
        Longitude of the center of the region.
    mol : str
        Molecule name.
    """

    reg = regions.RectangleSkyRegion(center=SkyCoord(l, (B_MIN+B_MAX)/2., frame='galactic'), width=1*u.arcmin, height=0.5*u.deg)
    subcube = cube.subcube_from_regions([reg])
    pv_mean = subcube.mean(axis=2)
    pv_mean.save(f'{save_path}/{mol}_pv_l{round(l.value, 3)}_mean.fits')
    pv_max = subcube.max(axis=2)
    pv_max.save(f'{save_path}/{mol}_pv_l{round(l.value, 3)}_max.fits')

def make_pv_mol(cube_fn):
    """
    Make PV diagrams for a given molecule.

    Parameters
    ----------
    cube_fn : str
        Path to the fits file of the cube.
    """
    
    cube = SpectralCube.read(cube_fn)
    mol = cube_fn.split('/')[-1].split('_')[0]

    list_b = make_position_list(B_MIN, B_MAX, 0.5*u.arcmin, u.arcmin)
    list_l = make_position_list(L_MIN, L_MAX, 0.5*u.arcmin, u.arcmin)

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