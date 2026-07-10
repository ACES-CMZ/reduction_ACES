"""
Summary figures for the ACES maser crossmatch (see crossmatch_masers.py).

Produces:

1. ``aces_maser_crossmatch_lb.png`` -- l/b scatter of all ACES compact
   continuum sources (low opacity) with maser-matched sources highlighted
   (H2O, CH3OH class I/II, SiO), within the match radius.
2. ``aces_continuum_masers_overlay.png`` -- the ACES 12m continuum mosaic with
   all H2O, CH3OH (class I/II) and SiO masers in the field overlaid.

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

# Maser types are kept distinct: H2O, CH3OH class I (36 GHz, collisionally
# pumped, trace shocks/outflows), CH3OH class II (6.7 GHz, radiatively pumped,
# trace HMYSOs), and SiO (evolved-star / stellar masers).  Each gets its own
# colour + marker in every figure.
MASER_TYPES = [
    ("H2O", "n_H2O_masers", "H2O", "#1f77b4", "o"),
    ("CH3OH class I (36 GHz)", "n_CH3OH_classI_masers", "CH3OH_classI", "#ff7f0e", "^"),
    ("CH3OH class II (6.7 GHz)", "n_CH3OH_classII_masers", "CH3OH_classII", "#2ca02c", "s"),
    ("SiO (stellar)", "n_SiO_masers", "SiO", "#d62728", "D"),
]


def _wrap_lon(lon):
    lon = np.asarray(lon, dtype=float)
    return np.where(lon > 180, lon - 360, lon)


def plot_lb_scatter(aces, outfile, radius_arcsec=DEFAULT_RADIUS_ARCSEC):
    """l/b scatter: all sources faint, maser-matched sources highlighted by type."""
    lon = _wrap_lon(aces["GLON_peak"])
    lat = np.asarray(aces["GLAT_peak"], dtype=float)

    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.scatter(lon, lat, s=14, c="0.6", alpha=0.25, edgecolors="none",
               label=f"All ACES sources (N={len(aces)})", zorder=1)
    for label, col, _mt, color, marker in MASER_TYPES:
        if col not in aces.colnames:
            continue
        sel = np.asarray(aces[col], dtype=int) > 0
        ax.scatter(lon[sel], lat[sel], s=60, c=color, alpha=0.9, marker=marker,
                   edgecolors="k", linewidths=0.4,
                   label=f"{label} (N={int(sel.sum())})", zorder=3)

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
    """ACES continuum mosaic with masers overlaid, by type (H2O, CH3OH I/II, SiO)."""
    with fits.open(continuum_fits) as hdul:
        data = np.squeeze(hdul[0].data)
        header = hdul[0].header
    ww = WCS(header).celestial

    maser_coords = SkyCoord(np.asarray(masers["ra_deg"]) * u.deg,
                            np.asarray(masers["dec_deg"]) * u.deg)
    maser_type = np.array([str(t) for t in masers["maser_type"]])

    # pixel coordinates; keep only masers that land within the mosaic
    px, py = ww.world_to_pixel(maser_coords.galactic)
    ny, nx = data.shape
    in_field = (px >= 0) & (px < nx) & (py >= 0) & (py < ny)

    fig = plt.figure(figsize=(15, 5.5))
    ax = fig.add_subplot(projection=ww)
    norm = simple_norm(data, stretch="asinh", min_percent=1, max_percent=99.7)
    ax.imshow(data, origin="lower", cmap="gray_r", norm=norm)

    for label, _col, mt, color, marker in MASER_TYPES:
        sel = in_field & (maser_type == mt)
        ax.scatter(px[sel], py[sel], s=45, facecolors="none", edgecolors=color,
                   linewidths=1.2, marker=marker,
                   label=f"{label} in field (N={int(sel.sum())})")

    ax.coords[0].set_axislabel("Galactic Longitude")
    ax.coords[1].set_axislabel("Galactic Latitude")
    ax.set_title("ACES 12m continuum with H2O, CH3OH (class I / II) and SiO masers overlaid")
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
