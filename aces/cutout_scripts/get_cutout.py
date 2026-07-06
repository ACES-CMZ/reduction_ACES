#!/usr/bin/env python3
"""Download an ALMA SODA cutout from Galactic coordinates."""

import argparse
from pathlib import Path


SODA_URL = "https://almascience.eso.org/soda/sync"


def galactic_to_icrs(glon, glat):
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    coord = SkyCoord(l=glon * u.deg, b=glat * u.deg, frame="galactic")
    return coord.icrs.ra.deg, coord.icrs.dec.deg


def output_name(filename, glon, glat, radius, vmin=None, vmax=None):
    stem = Path(filename).name.removesuffix(".fits")
    if vmin is not None and vmax is not None:
        return (
            f"{stem}_l{glon:g}_b{glat:g}_r{radius:g}arcsec_"
            f"v{vmin:g}_{vmax:g}kms_cutout.fits"
        )
    return f"{stem}_l{glon:g}_b{glat:g}_r{radius:g}arcsec_cutout.fits"


def fits_value(card):
    value = card[10:30].split("/", 1)[0].strip()
    if value.startswith("'"):
        return value.strip("' ")
    try:
        return float(value)
    except ValueError:
        return value


def read_primary_header(response):
    header = b""

    for chunk in response.iter_content(chunk_size=2880):
        header += chunk
        try:
            text = header.decode("ascii")
        except UnicodeDecodeError as exc:
            raise RuntimeError("SODA response did not start with a FITS header") from exc

        cards = [text[i:i + 80] for i in range(0, len(text), 80)]
        if any(card.startswith("END") for card in cards):
            return {
                card[:8].strip(): fits_value(card)
                for card in cards
                if "=" in card[:10]
            }

    raise RuntimeError("Could not find END card in FITS header")


def fetch_primary_header(requests_module, params):
    response = requests_module.get(SODA_URL, params=params, stream=True, timeout=120)
    content_type = response.headers.get("content-type", "")

    if "fits" not in content_type.lower():
        print(response.url)
        print(response.status_code)
        print(content_type)
        print(response.text[:2000])
        response.raise_for_status()
        raise RuntimeError("SODA header response did not look like a FITS file")

    try:
        return read_primary_header(response)
    finally:
        response.close()


def velocity_to_band(restfrq_hz, vmin, vmax):
    from astropy import constants as const
    from astropy import units as u

    rest_frequency = restfrq_hz * u.Hz
    wavelengths = []
    for velocity_kms in (vmin, vmax):
        velocity = velocity_kms * u.km / u.s
        frequency = rest_frequency * (1.0 - (velocity / const.c).decompose())
        wavelength = (const.c / frequency).to(u.m)
        wavelengths.append(wavelength.value)
    return sorted(wavelengths)


def parse_extra_args(parser, extra_args):
    if len(extra_args) == 0:
        return None, None, None
    if len(extra_args) == 1:
        return None, None, extra_args[0]
    if len(extra_args) == 2:
        try:
            return float(extra_args[0]), float(extra_args[1]), None
        except ValueError:
            parser.error("vmin and vmax must be numbers in km/s")
    if len(extra_args) == 3:
        try:
            return float(extra_args[0]), float(extra_args[1]), extra_args[2]
        except ValueError:
            parser.error("vmin and vmax must be numbers in km/s")

    parser.error("expected: filename glon glat radius [vmin vmax] [outfile]")


def get_cutout(filename, glon, glat, radius, vmin=None, vmax=None, outfile=None):
    from astropy import units as u
    import requests

    ra, dec = galactic_to_icrs(glon, glat)
    radius_deg = (radius * u.arcsec).to_value(u.deg)

    params = {
        "REQUEST": "queryData",
        "ID": filename,
        "POS": f"CIRCLE {ra:.8f} {dec:.8f} {radius_deg:.8f}",
    }

    if vmin is not None and vmax is not None:
        header = fetch_primary_header(requests, params)
        restfrq = header.get("RESTFRQ", header.get("RESTFREQ"))
        if restfrq is None:
            raise RuntimeError("Could not find RESTFRQ or RESTFREQ in FITS header")

        band_min, band_max = velocity_to_band(float(restfrq), vmin, vmax)
        params["BAND"] = f"{band_min:.8e} {band_max:.8e}"

    response = requests.get(SODA_URL, params=params, timeout=120)
    content_type = response.headers.get("content-type", "")

    print(response.url)
    print(response.status_code)
    print(content_type)

    if "fits" not in content_type.lower():
        print(response.text[:2000])
        response.raise_for_status()
        raise RuntimeError("SODA response did not look like a FITS file")

    if outfile is None:
        outfile = output_name(filename, glon, glat, radius, vmin, vmax)

    Path(outfile).write_bytes(response.content)
    print(f"Wrote {outfile}")


def main():
    parser = argparse.ArgumentParser(
        description="Download an ALMA SODA cutout using Galactic coordinates."
    )
    parser.add_argument("filename", help="SODA dataset ID or FITS filename")
    parser.add_argument("glon", type=float, help="Galactic longitude in degrees")
    parser.add_argument("glat", type=float, help="Galactic latitude in degrees")
    parser.add_argument("radius", type=float, help="Cutout radius in arcsec")
    parser.add_argument(
        "extra_args",
        nargs="*",
        help="Optional [vmin vmax] in km/s followed by optional output FITS filename",
    )
    args = parser.parse_args()
    vmin, vmax, outfile = parse_extra_args(parser, args.extra_args)

    get_cutout(args.filename, args.glon, args.glat, args.radius, vmin, vmax, outfile)


if __name__ == "__main__":
    main()
