# ALMA SODA Cutout Script

This directory contains a small helper script, `get_cutout.py`, for downloading
ALMA Science Archive SODA cutouts from Galactic coordinates.

The script accepts a SODA dataset ID, a Galactic longitude and latitude, and a
spatial radius. For spectral cubes, it can also accept a radio velocity range
and translate that range into the wavelength `BAND` parameter expected by SODA.

## Requirements

Use a Python environment with:

- `requests`
- `astropy`
- `matplotlib`
- `numpy`

Install them with:

```bash
python -m pip install -r requirements.txt
```

Check the command interface with:

```bash
python get_cutout.py --help
```

If that Python environment does not have the required packages, activate an
environment that includes the packages in `requirements.txt`. On this machine,
the Anaconda Python at `/opt/anaconda3/bin/python` has these dependencies
available.

The `cutout_examples.ipynb` notebook uses the same requirements and shows image
and cube examples using this script as an import.

## Basic Usage

General syntax:

```bash
python get_cutout.py filename glon glat radius [vmin vmax] [outfile]
```

Arguments:

- `filename`: ALMA SODA dataset ID or archive FITS filename.
- `glon`: Galactic longitude in degrees.
- `glat`: Galactic latitude in degrees.
- `radius`: Spatial cutout radius in arcsec.
- `vmin`: Optional minimum radio velocity in `km/s`.
- `vmax`: Optional maximum radio velocity in `km/s`.
- `outfile`: Optional output FITS filename.

You can find the available ACES group-level FITS filenames on the ALMA Science
Archive page:

https://almascience.eso.org/alma-data/lp/aces-group-level-data

This directory also includes `known_aces_files.dat`, a local list of known ACES
FITS filenames that can be compared against when choosing or checking a
`filename` value.

The supported command forms are:

```bash
python get_cutout.py filename glon glat radius
python get_cutout.py filename glon glat radius output.fits
python get_cutout.py filename glon glat radius vmin vmax
python get_cutout.py filename glon glat radius vmin vmax output.fits
```

If `outfile` is not provided, the script writes a generated filename based on
the dataset ID, position, radius, and optional velocity range.

## Examples

### Spatial Cutout From A Continuum Mosaic

This downloads a 20 arcsec radius cutout around Galactic coordinates
`l=0.0`, `b=0.0` from the 93.7 GHz continuum mosaic and writes the result to
`sgra_test.fits`.

```bash
python get_cutout.py \
  group.uid___A001_X1590_X30a9.lp_slongmore.cmz_mosaic.icrs.12m.cont.93.7GHz_bw4.7GHz.pbcor.fits \
  0.0 0.0 20 \
  sgra_test.fits
```

The same command without an explicit output filename will write an automatic
name:

```bash
python get_cutout.py \
  group.uid___A001_X1590_X30a9.lp_slongmore.cmz_mosaic.icrs.12m.cont.93.7GHz_bw4.7GHz.pbcor.fits \
  0.0 0.0 20
```

### Velocity Cutout From The H40a Cube

This downloads a 20 arcsec radius cutout around `l=0.0`, `b=0.0`, limited to
radio velocities from `-100` to `100 km/s`, and writes the result to
`h40a_test.fits`.

```bash
python get_cutout.py \
  group.uid___A001_X1590_X30a9.lp_slongmore.cmz_mosaic.icrs.12m7mTP.H40a.cube.pbcor.fits \
  0.0 0.0 20 \
  -100 100 \
  h40a_test.fits
```

A successful velocity cutout request prints a URL containing both the spatial
cutout and the spectral `BAND` interval, for example:

```text
POS=CIRCLE+266.40498829+-28.93617776+0.00555556&BAND=3.02649529e-03+3.02851502e-03
```

## What The Script Does

For every request, the script:

1. Converts `glon` and `glat` from Galactic coordinates to ICRS RA and Dec.
2. Converts `radius` from arcsec to degrees.
3. Builds a SODA `POS=CIRCLE ra dec radius_deg` request.
4. Downloads the cutout from `https://almascience.eso.org/soda/sync`.
5. Writes the returned FITS content to the requested output file.

For velocity cutouts, the script does one extra preliminary request:

1. It streams the beginning of the FITS response.
2. It reads only the primary FITS header.
3. It extracts `RESTFRQ`, falling back to `RESTFREQ`.
4. It converts the requested radio velocities in `km/s` to wavelength limits in
   meters.
5. It adds those limits as the SODA `BAND` parameter.

The radio velocity conversion is:

```text
nu = restfrq_hz * (1 - velocity_kms * 1000 / c)
lambda_m = c / nu
```

The two wavelengths are sorted before being passed to SODA as `BAND`.

The script does not transform velocity frames. It assumes the requested
velocities are in the same practical frame as the archive cube header and rest
frequency.

## Output Printed By The Script

The script prints:

- The final SODA URL used for the download.
- The HTTP status code.
- The response content type.
- The output FITS filename after writing.

Example:

```text
https://almascience.eso.org/soda/sync?REQUEST=queryData&ID=...&POS=...&BAND=...
200
application/fits
Wrote h40a_test.fits
```

If the response is not FITS, the script prints the first part of the service
response text before raising an error. This is useful for archive-side errors,
bad dataset IDs, unsupported cutout requests, or authentication/service issues.

## Troubleshooting

### `ModuleNotFoundError: No module named 'requests'`

Use a Python environment with `requests` and `astropy`:

```bash
python get_cutout.py --help
```

If `python` points to an environment without those packages, activate your
science Python environment or use `/opt/anaconda3/bin/python get_cutout.py`.

### The Response Is Not FITS

Check the printed URL, status code, content type, and response text. Common
causes are:

- The `filename` is not a valid SODA dataset ID.
- The requested spatial region is outside the product footprint.
- The requested velocity range is outside the product spectral coverage.
- The archive service returned an error page or diagnostic message.

### Velocity Cutouts Fail With Missing Rest Frequency

Velocity cutouts require `RESTFRQ` or `RESTFREQ` in the FITS header. If neither
keyword is present, the script cannot convert velocity limits into a SODA
`BAND` wavelength interval.

In that case, either use a spatial-only cutout or update the script to accept an
explicit rest frequency for that product.

### Network Or DNS Errors

The script needs network access to `almascience.eso.org`. If a sandboxed
environment reports DNS or connection errors, rerun the command in a normal
network-enabled shell.

### Checking The Command Interface

Run:

```bash
python get_cutout.py --help
```

Expected output:

```text
usage: get_cutout.py [-h] filename glon glat radius [extra_args ...]

Download an ALMA SODA cutout using Galactic coordinates.

positional arguments:
  filename    SODA dataset ID or FITS filename
  glon        Galactic longitude in degrees
  glat        Galactic latitude in degrees
  radius      Cutout radius in arcsec
  extra_args  Optional [vmin vmax] in km/s followed by optional output FITS
              filename

options:
  -h, --help  show this help message and exit
```
