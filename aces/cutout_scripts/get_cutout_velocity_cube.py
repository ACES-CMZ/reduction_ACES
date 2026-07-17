#!/usr/bin/env python3
"""Download a spectral cutout and save it as a velocity cube in Kelvin."""

import argparse

from astropy import constants as const
from astropy import units as u

import get_cutout


def positive_number(value):
    """Convert a command-line value to a positive number."""
    number = float(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return number


def frequency_limits(rfreq, vrange):
    """Return the frequencies for -vrange and +vrange around a spectral line.

    Parameters
    ----------
    rfreq : float
        Rest frequency in GHz.
    vrange : float
        Positive velocity half-range in km/s.
    """
    rest_frequency = rfreq * u.GHz
    velocities = (-vrange, vrange) * u.km / u.s

    # Radio velocity convention: frequency = rest_frequency * (1 - velocity/c)
    frequencies = rest_frequency * (1 - velocities / const.c)
    return sorted(frequencies.to_value(u.GHz))


def get_cutout_velocity_cube(filename, glon, glat, radius, rfreq, vrange, outfile):
    """Download and process one spectral-line cube."""
    from spectral_cube import SpectralCube

    fmin, fmax = frequency_limits(rfreq, vrange)

    print(f"Downloading frequencies from {fmin:.6f} to {fmax:.6f} GHz")
    get_cutout.get_cutout(
        filename,
        glon,
        glat,
        radius,
        fmin=fmin,
        fmax=fmax,
        outfile=outfile,
    )

    print("Converting the downloaded cube to radio velocity and Kelvin")
    cube = SpectralCube.read(outfile)
    cube.allow_huge_operations = True
    cube = cube.with_spectral_unit(
        u.km / u.s,
        velocity_convention="radio",
        rest_value=rfreq * u.GHz,
    )
    cube.allow_huge_operations = True
    cube = cube.to(u.K)

    # Replace the SODA download, so only the finished cube remains.
    cube.write(outfile, overwrite=True)
    print(f"Finished: {outfile}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download a SODA cutout around a spectral line and save it with "
            "a radio-velocity axis in km/s and brightness units of K."
        )
    )
    parser.add_argument("filename", help="SODA dataset ID or archive FITS filename")
    parser.add_argument("glon", type=float, help="Galactic longitude in degrees")
    parser.add_argument("glat", type=float, help="Galactic latitude in degrees")
    parser.add_argument(
        "radius", type=positive_number, help="Cutout radius in arcsec"
    )
    parser.add_argument(
        "rfreq", type=positive_number, help="Line rest frequency in GHz"
    )
    parser.add_argument(
        "vrange",
        type=positive_number,
        help="Velocity half-range in km/s (for example, 100 means -100 to +100)",
    )
    parser.add_argument("outfile", help="Output FITS filename")
    args = parser.parse_args()

    get_cutout_velocity_cube(
        args.filename,
        args.glon,
        args.glat,
        args.radius,
        args.rfreq,
        args.vrange,
        args.outfile,
    )


if __name__ == "__main__":
    main()
