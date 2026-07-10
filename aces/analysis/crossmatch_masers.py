"""
Crossmatch the ACES compact continuum source catalog against maser catalogs.

Included maser catalogs (bundled under ``aces/data/masers/``; regenerate with
prepare_maser_catalogs.py):

* **H2O** -- Walsh et al. (2014) ATCA water masers (all-plane).
* **H2O** -- Lu et al. (2019) CMZ water masers.

* **CH3OH** -- Cotton & Yusef-Zadeh (2016) VLA 36 GHz class-I masers, the dense
  CMZ methanol catalog (2240 spots).
* **CH3OH** -- GLOSTAR 6.7 GHz class-II masers (Nguyen et al. 2022, all-plane).

Additionally read from disk at runtime when present (not redistributed):

* **H2O** -- SWAG 22 GHz water masers (Ward, Ott & Meier, in prep.; CMZ;
  ``SWAG_DIR``).
* **SiO** -- ACES evolved-star (stellar) SiO masers detected in the ACES data
  itself (``SIO_FILE``), collapsed to one entry per source.

For every ACES source we record, within a configurable match radius:

* ``maser_match``            -- True if any maser (any type) is within radius
* ``n_H2O_masers``, ``n_CH3OH_classI_masers``, ``n_CH3OH_classII_masers``
  -- number of masers within radius, counted separately by maser type.
  CH3OH class I (36 GHz) and class II (6.7 GHz) are kept distinct throughout.
* ``maser_nearest_*``        -- properties of the single nearest maser
  (separation, catalog, species, class, type, name, V_lsr, flux) when within radius

Outputs ``aces_compact_catalog_v0_withMasers.{fits,ecsv}`` and a per-match
table ``aces_maser_matches.ecsv`` under ``{basepath}/tables/``.

Run: ``aces_crossmatch_masers`` (or ``python -m aces.analysis.crossmatch_masers``).
"""
import argparse
import os
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic, search_around_sky
from astropy.table import Table, vstack

from aces import conf

basepath = Path(conf.basepath)

# ACES compact continuum source catalog (Galactic coords, key column "index")
DEFAULT_CATALOG_PATH = basepath / "tables" / "aces_compact_catalog_v0_withExtras.fits"
OUTPUT_DIR = basepath / "tables"

# Maser catalogs bundled with the package
MASER_DATA_DIR = Path(__file__).parent.parent / "data" / "masers"

# SWAG water masers (Ward, Ott & Meier, in prep.) are NOT redistributed; they are
# read from disk at runtime if the machine-readable tables are present.
SWAG_DIR = Path("/orange/adamginsburg/cmz/swag/dylanward_watermasers/Tables")

# ACES SiO masers -- evolved-star (stellar) SiO masers detected in the ACES data
# itself.  These are an ACES data product, not a redistributed external catalog,
# so like SWAG they are read from disk at runtime when present.  The
# position_velocity table carries one row per fitted spectral peak; it is
# collapsed to one row per source in _load_sio_masers.
SIO_FILE = basepath / "upload" / "SiO_Masers" / "ACES_SiO_maser_position_velocity.fits"

# Default association radius between a continuum source and a maser
DEFAULT_RADIUS_ARCSEC = 2.0

# Bundled, pre-normalized maser catalogs (regenerate with
# prepare_maser_catalogs.py).  Each file has uniform columns:
# ra_deg, dec_deg, vlsr, flux, flux_unit, name, species, catalog, reference.
MASER_CATALOG_FILES = [
    "h2o_walsh2014.fits",        # H2O, ATCA all-plane (Walsh et al. 2014)
    "h2o_lu2019_cmz.fits",       # H2O, CMZ (Lu et al. 2019)
    "ch3oh_cotton2016_cmz.fits",  # CH3OH 36 GHz class I, CMZ (Cotton & Yusef-Zadeh 2016)
    "ch3oh_glostar.fits",        # CH3OH 6.7 GHz class II, all-plane (GLOSTAR)
]


def load_aces_catalog(catalog_path=DEFAULT_CATALOG_PATH):
    """Load the ACES catalog and attach equatorial coordinates.

    The catalog stores Galactic ``GLON_peak``/``GLAT_peak`` (deg); we add
    ``RA_J2000``/``Dec_J2000`` and ensure an ``index`` join key exists.
    """
    cat = Table.read(catalog_path)
    gc = SkyCoord(l=cat["GLON_peak"], b=cat["GLAT_peak"], frame=Galactic, unit="deg")
    icrs = gc.icrs
    cat["RA_J2000"] = icrs.ra.deg
    cat["Dec_J2000"] = icrs.dec.deg
    if "index" not in cat.colnames:
        cat["index"] = np.arange(len(cat))
    return cat


def _load_swag_masers(swag_dir=SWAG_DIR):
    """Read the SWAG water maser catalog (Ward, Ott & Meier, in prep.) from its
    machine-readable tables on disk and normalize it to the bundled schema.

    Returns None if the tables are not present (so the crossmatch runs without
    it); the catalog is not redistributed with the package.
    """
    swag_dir = Path(swag_dir)
    pos_path = swag_dir / "MR_astrodendro_ironclad_updated.txt"
    fit_path = swag_dir / "MR_fits_updated.txt"
    if not (pos_path.exists() and fit_path.exists()):
        print(f"SWAG water masers not found at {swag_dir}; skipping (not bundled).")
        return None

    pos = Table.read(pos_path, format="ascii.cds")
    fts = Table.read(fit_path, format="ascii.cds")
    # per source, keep the strongest Gaussian component (blank ID rows are
    # continuation components of the previous source)
    fmask = fts["ID"].mask if hasattr(fts["ID"], "mask") else np.zeros(len(fts), bool)
    amp = np.asarray(fts["Amp"], float)
    cen = np.asarray(fts["Centroid"], float)
    best = {}
    last = -1
    for i in range(len(fts)):
        if not fmask[i]:
            last = int(fts["ID"][i])
        if last not in best or amp[i] > best[last][0]:
            best[last] = (amp[i], cen[i])
    pid = np.asarray(pos["ID"], int)
    coords = SkyCoord(l=np.asarray(pos["GLON"], float) * u.deg,
                      b=np.asarray(pos["GLAT"], float) * u.deg, frame=Galactic).icrs

    out = Table()
    out["ra_deg"] = coords.ra.deg
    out["dec_deg"] = coords.dec.deg
    out["vlsr"] = np.array([best.get(s, (np.nan, np.nan))[1] for s in pid])
    out["flux"] = np.array([best.get(s, (np.nan, np.nan))[0] for s in pid])
    out["flux_unit"] = "K"
    out["name"] = np.array([f"SWAG-{s}" for s in pid])
    out["species"] = "H2O"
    out["maser_class"] = ""
    out["maser_type"] = "H2O"
    out["catalog"] = "H2O_SWAG_Ward"
    out["reference"] = "Ward, Ott & Meier (SWAG 22 GHz water masers, CMZ; in prep.)"
    return out


def _load_sio_masers(sio_file=SIO_FILE):
    """Read the ACES SiO (evolved-star / stellar) maser catalog from disk and
    normalize it to the bundled schema.

    Returns None if the table is not present (so the crossmatch runs without it);
    it is an ACES data product and is not redistributed with the package.  The
    position_velocity table has one row per fitted spectral peak, so collapse it
    to one row per source (the position is identical across a source's peaks;
    keep the first peak's velocity as the representative V_lsr).
    """
    sio_file = Path(sio_file)
    if not sio_file.exists():
        print(f"ACES SiO masers not found at {sio_file}; skipping (not bundled).")
        return None

    t = Table.read(sio_file)
    t = t.group_by("number")
    src = t[t.groups.indices[:-1]]  # first row (peak) of each source group
    coords = SkyCoord(l=np.asarray(src["GLON"], float) * u.deg,
                      b=np.asarray(src["GLAT"], float) * u.deg, frame=Galactic).icrs

    out = Table()
    out["ra_deg"] = coords.ra.deg
    out["dec_deg"] = coords.dec.deg
    out["vlsr"] = np.asarray(src["velocity"], float)
    out["flux"] = np.full(len(src), np.nan)
    out["flux_unit"] = ""
    out["name"] = np.array([f"ACES-SiO-{int(n)}" for n in src["number"]])
    out["species"] = "SiO"
    out["maser_class"] = "stellar"
    out["maser_type"] = "SiO"
    out["catalog"] = "SiO_ACES"
    out["reference"] = "ACES SiO maser survey (evolved-star SiO masers, CMZ; this work)"
    return out


def load_masers(data_dir=MASER_DATA_DIR, files=MASER_CATALOG_FILES, include_swag=True,
                include_sio=True):
    """Load and concatenate the bundled maser catalogs, plus the SWAG water
    masers and ACES SiO masers read from disk at runtime when available (neither
    is bundled)."""
    tables = []
    for fn in files:
        path = Path(data_dir) / fn
        if not path.exists():
            print(f"WARNING: maser catalog {path} not found; skipping")
            continue
        tbl = Table.read(path)
        species = str(tbl["species"][0]) if len(tbl) else "?"
        print(f"Loaded {len(tbl)} masers from {fn} ({species})")
        tables.append(tbl)
    if include_swag:
        swag = _load_swag_masers()
        if swag is not None:
            print(f"Loaded {len(swag)} masers from SWAG water masers (disk, H2O)")
            tables.append(swag)
    if include_sio:
        sio = _load_sio_masers()
        if sio is not None:
            print(f"Loaded {len(sio)} masers from ACES SiO masers (disk, SiO stellar)")
            tables.append(sio)
    if not tables:
        raise RuntimeError("No maser catalogs could be loaded")
    return vstack(tables, metadata_conflicts="silent")


def crossmatch_masers(aces, masers, radius_arcsec=DEFAULT_RADIUS_ARCSEC):
    """Augment the ACES catalog with maser-association columns and return
    (augmented_catalog, per-match table)."""
    radius = radius_arcsec * u.arcsec
    aces_coords = SkyCoord(aces["RA_J2000"] * u.deg, aces["Dec_J2000"] * u.deg)
    maser_coords = SkyCoord(masers["ra_deg"] * u.deg, masers["dec_deg"] * u.deg)

    nsrc = len(aces)
    # count/report by maser_type, which distinguishes CH3OH class I vs II
    # (H2O, CH3OH_classI, CH3OH_classII).
    maser_type = np.array([str(t) for t in masers["maser_type"]])
    type_list = sorted(set(maser_type))

    # nearest maser overall
    idx_near, sep_near, _ = aces_coords.match_to_catalog_sky(maser_coords)
    within = sep_near < radius

    aces = aces.copy()
    aces["maser_match"] = within
    aces["maser_nearest_sep_arcsec"] = np.where(within, sep_near.arcsec, np.nan)
    for col, key, fill in [("maser_nearest_catalog", "catalog", ""),
                           ("maser_nearest_species", "species", ""),
                           ("maser_nearest_class", "maser_class", ""),
                           ("maser_nearest_type", "maser_type", ""),
                           ("maser_nearest_name", "name", ""),
                           ("maser_nearest_reference", "reference", "")]:
        vals = np.array([str(masers[key][i]) for i in idx_near])
        aces[col] = np.where(within, vals, fill)
    for col, key in [("maser_nearest_vlsr", "vlsr"), ("maser_nearest_flux", "flux")]:
        vals = np.asarray(masers[key])[idx_near]
        aces[col] = np.where(within, vals, np.nan)

    # counts within radius, per maser_type (search_around_sky finds all pairs).
    # Function form: ia indexes aces_coords, im indexes maser_coords.
    ia, im, sep2d, _ = search_around_sky(aces_coords, maser_coords, radius)
    im_type = maser_type[im]
    for mt in type_list:
        counts = np.zeros(nsrc, dtype=int)
        sel = im_type == mt
        if sel.any():
            np.add.at(counts, ia[sel], 1)
        aces[f"n_{mt}_masers"] = counts

    # per-match table (one row per ACES-maser pair within radius)
    if len(ia):
        match = Table()
        match["index"] = np.asarray(aces["index"])[ia]
        match["sep_arcsec"] = sep2d.arcsec
        for key in ("catalog", "species", "maser_class", "maser_type", "name", "vlsr", "flux", "reference"):
            match[f"maser_{key}" if not key.startswith("maser_") else key] = np.asarray(masers[key])[im]
        match["maser_ra_deg"] = np.asarray(masers["ra_deg"])[im]
        match["maser_dec_deg"] = np.asarray(masers["dec_deg"])[im]
        match.sort(["index", "sep_arcsec"])
    else:
        match = Table(names=["index", "sep_arcsec", "maser_catalog", "maser_species",
                             "maser_class", "maser_type", "maser_name", "maser_vlsr",
                             "maser_flux", "maser_reference", "maser_ra_deg", "maser_dec_deg"])

    return aces, match


def main(argv=None):
    parser = argparse.ArgumentParser(description="Crossmatch ACES sources against H2O/CH3OH maser catalogs")
    parser.add_argument("--catalog", default=str(DEFAULT_CATALOG_PATH),
                        help="Path to the ACES base catalog (FITS)")
    parser.add_argument("--radius", type=float,
                        default=float(os.environ.get("MASER_MATCH_RADIUS", DEFAULT_RADIUS_ARCSEC)),
                        help="Match radius in arcsec (default 2)")
    parser.add_argument("--outdir", default=str(OUTPUT_DIR),
                        help="Output directory")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading ACES catalog from {args.catalog}")
    aces = load_aces_catalog(args.catalog)
    print(f"  {len(aces)} sources")

    masers = load_masers()
    print(f"Total masers loaded: {len(masers)}")

    aces, match = crossmatch_masers(aces, masers, radius_arcsec=args.radius)

    nmatched = int(np.count_nonzero(aces["maser_match"]))
    print(f"\n{nmatched}/{len(aces)} ACES sources have a maser within {args.radius}\"")
    for mt in sorted(set(str(t) for t in masers["maser_type"])):
        col = f"n_{mt}_masers"
        if col in aces.colnames:
            print(f"  sources with >=1 {mt} maser: {int(np.count_nonzero(aces[col]))}")
    print(f"  total ACES-maser pairs within radius: {len(match)}")

    out_fits = outdir / "aces_compact_catalog_v0_withMasers.fits"
    out_ecsv = outdir / "aces_compact_catalog_v0_withMasers.ecsv"
    match_ecsv = outdir / "aces_maser_matches.ecsv"
    aces.write(out_fits, overwrite=True)
    aces.write(out_ecsv, format="ascii.ecsv", overwrite=True)
    match.write(match_ecsv, format="ascii.ecsv", overwrite=True)
    print(f"\nWrote {out_fits}\n      {out_ecsv}\n      {match_ecsv}")


if __name__ == "__main__":
    main()
