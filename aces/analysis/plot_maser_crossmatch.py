"""
Summary figures for the ACES maser crossmatch (see crossmatch_masers.py).

Produces:

1. ``aces_maser_crossmatch_lb.png`` -- l/b scatter of all ACES compact
   continuum sources (low opacity) with maser-matched sources highlighted
   (H2O, CH3OH, both), within the match radius.
2. ``aces_continuum_masers_overlay.png`` -- the ACES 12m continuum mosaic with
   all H2O and CH3OH masers in the field overlaid.

Run: ``aces_plot_maser_crossmatch`` (or
``python -m aces.analysis.plot_maser_crossmatch``).
"""
import argparse
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import simple_norm

from aces import conf
from aces.analysis.crossmatch_masers import (
    load_aces_catalog, load_masers, crossmatch_masers, DEFAULT_RADIUS_ARCSEC,
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

basepath = Path(conf.basepath)

DEFAULT_CONTINUUM = (basepath / "mosaics" / "continuum" /
                     "12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits")
DEFAULT_OUTDIR = basepath / "figures"

H2O_COLOR = "#1f77b4"
CH3OH_COLOR = "#d62728"
BOTH_COLOR = "#9467bd"


def _wrap_lon(lon):
    lon = np.asarray(lon, dtype=float)
    return np.where(lon > 180, lon - 360, lon)


def plot_lb_scatter(aces, outfile, radius_arcsec=DEFAULT_RADIUS_ARCSEC):
    """l/b scatter: all sources faint, maser-matched sources highlighted."""
    lon = _wrap_lon(aces["GLON_peak"])
    lat = np.asarray(aces["GLAT_peak"], dtype=float)
    h2o = np.asarray(aces["n_H2O_masers"], dtype=int) > 0
    ch3oh = np.asarray(aces["n_CH3OH_masers"], dtype=int) > 0
    both = h2o & ch3oh

    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.scatter(lon, lat, s=14, c="0.6", alpha=0.25, edgecolors="none",
               label=f"All ACES sources (N={len(aces)})", zorder=1)
    ax.scatter(lon[h2o & ~both], lat[h2o & ~both], s=55, c=H2O_COLOR, alpha=0.95,
               edgecolors="k", linewidths=0.4,
               label=f"H2O maser (N={int((h2o & ~both).sum())})", zorder=3)
    ax.scatter(lon[ch3oh & ~both], lat[ch3oh & ~both], s=55, c=CH3OH_COLOR, alpha=0.95,
               edgecolors="k", linewidths=0.4,
               label=f"CH3OH maser (N={int((ch3oh & ~both).sum())})", zorder=3)
    ax.scatter(lon[both], lat[both], s=75, marker="D", c=BOTH_COLOR, alpha=0.98,
               edgecolors="k", linewidths=0.5,
               label=f"H2O + CH3OH (N={int(both.sum())})", zorder=4)

    ax.set_xlabel("Galactic Longitude [deg]")
    ax.set_ylabel("Galactic Latitude [deg]")
    ax.set_title(f"ACES compact continuum sources: maser crossmatch ({radius_arcsec:g} arcsec radius)")
    ax.invert_xaxis()
    ax.legend(loc="upper right", framealpha=0.9, fontsize=9)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outfile}")


def plot_continuum_overlay(masers, outfile, continuum_fits=DEFAULT_CONTINUUM):
    """ACES continuum mosaic with H2O and CH3OH masers overlaid."""
    with fits.open(continuum_fits) as hdul:
        data = np.squeeze(hdul[0].data)
        header = hdul[0].header
    ww = WCS(header).celestial

    maser_coords = SkyCoord(np.asarray(masers["ra_deg"]) * u.deg,
                            np.asarray(masers["dec_deg"]) * u.deg)
    species = np.array([str(s) for s in masers["species"]])

    # pixel coordinates; keep only masers that land within the mosaic
    px, py = ww.world_to_pixel(maser_coords.galactic)
    ny, nx = data.shape
    in_field = (px >= 0) & (px < nx) & (py >= 0) & (py < ny)

    fig = plt.figure(figsize=(15, 5.5))
    ax = fig.add_subplot(projection=ww)
    norm = simple_norm(data, stretch="asinh", min_percent=1, max_percent=99.7)
    ax.imshow(data, origin="lower", cmap="gray_r", norm=norm)

    for label, color, marker, sp in [("H2O maser", H2O_COLOR, "o", "H2O"),
                                     ("CH3OH maser", CH3OH_COLOR, "^", "CH3OH")]:
        sel = in_field & (species == sp)
        ax.scatter(px[sel], py[sel], s=45, facecolors="none", edgecolors=color,
                   linewidths=1.2, marker=marker,
                   label=f"{label} in field (N={int(sel.sum())})")

    ax.coords[0].set_axislabel("Galactic Longitude")
    ax.coords[1].set_axislabel("Galactic Latitude")
    ax.set_title("ACES 12m continuum with H2O and CH3OH masers overlaid")
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    ax.legend(loc="upper right", framealpha=0.9, fontsize=10)
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outfile}")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Plot ACES maser crossmatch summary figures")
    parser.add_argument("--radius", type=float, default=DEFAULT_RADIUS_ARCSEC,
                        help="Match radius in arcsec (default 2)")
    parser.add_argument("--continuum", default=str(DEFAULT_CONTINUUM),
                        help="ACES continuum mosaic FITS for the overlay")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR), help="Output directory")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    aces = load_aces_catalog()
    masers = load_masers()
    aces, _ = crossmatch_masers(aces, masers, radius_arcsec=args.radius)

    plot_lb_scatter(aces, outdir / "aces_maser_crossmatch_lb.png", radius_arcsec=args.radius)
    plot_continuum_overlay(masers, outdir / "aces_continuum_masers_overlay.png",
                           continuum_fits=args.continuum)


if __name__ == "__main__":
    main()
