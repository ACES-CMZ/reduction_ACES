#!/usr/bin/env python
"""
Crossmatch the ACES compact source catalog with SIMBAD and multi-wavelength
catalogs from Vizier (2MASS, Spitzer/GLIMPSE, MeerKAT, VLA, X-ray) to
classify each source.

Uses astroquery's CDS XMatch service for efficient, vectorized catalog
crossmatching against Vizier catalogs, and SIMBAD TAP queries for bulk
object-type lookups.  Astroquery's built-in caching is used throughout
(no custom cache layer).

Classification logic
--------------------
For each ACES source the script assigns a best-guess ``object_type``:

* **dust**          – detected at 1 mm (CMZoom) *and* 1 mm–3 mm spectral
                      index α ≥ 2  →  thermal dust emission dominates.
* **evolved_star**  – bright 2MASS (K < 10) *or* GLIMPSE (4.5 µm < 9)
                      counterpart  →  likely luminous post-main-sequence
                      star (AGB / RSG / WR).
* **hii_region**    – radio counterpart (VLA or MeerKAT) without a bright
                      NIR stellar match  →  free-free from an H II region
                      powered by a massive star.
* **radio+nir**     – radio *and* bright NIR counterpart  →  massive star
                      with both ionised gas and a detectable photosphere.
* **xray_source**   – Chandra or Muno counterpart  →  unusual; flagged for
                      further investigation.
* **simbad:<otype>** – SIMBAD classification when no other rule fires.
* **unclassified**  – no multi-wavelength counterpart at all.  These are
                      the "new discoveries" of greatest interest.

Output
------
``<output_dir>/aces_crossmatched_catalog.fits``
    FITS table with one row per ACES source and one column per external
    catalog giving the matched source name / ID (or empty), plus the
    ``object_type`` column.

``<output_dir>/aces_crossmatched_catalog.ecsv``
    Same table in human-readable ECSV format.

Usage
-----
Run with the specified environment::

    /blue/adamginsburg/adamginsburg/miniconda3/envs/python313/bin/python \\
        crossmatch_catalog.py
"""

import numpy as np
from pathlib import Path

from astropy.table import Table, Column, join, hstack
from astropy.coordinates import SkyCoord, Galactic
from astropy import units as u

from astroquery.xmatch import XMatch
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
CATALOG_PATH = Path(
    "/orange/adamginsburg/ACES/tables/aces_compact_catalog_v0_withExtras.fits"
)
OUTPUT_DIR = Path("/orange/adamginsburg/ACES/tables")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Match radii (arcsec)  – chosen to reflect the angular resolution of each
# survey relative to the ~2″ ACES beam.
# ---------------------------------------------------------------------------
MATCH_RADIUS = {
    "2MASS":       3.0,  # 2″ 2MASS pixels; allow small offset
    "GLIMPSE":     3.0,  # ~2″ Spitzer/IRAC PSF
    "MeerKAT":     5.0,  # ~7″ MeerKAT beam; be generous
    "VLA_Lu2019":  3.0,  # ~2″ VLA CMZ data
    "Chandra_CSC": 3.0,  # sub-arcsec Chandra, but positional uncertainties
    "Muno_Xray":   3.0,  # Muno+ 2009 positional accuracy ~1″
}

# ---------------------------------------------------------------------------
# Vizier catalog identifiers
# ---------------------------------------------------------------------------
VIZIER_CATALOGS = {
    # 2MASS All-Sky Point Source Catalog (Cutri+ 2003)
    "2MASS":       "II/246/out",
    # Spitzer/GLIMPSE Source Catalog (I + II + 3D; IPAC 2008)
    "GLIMPSE":     "II/293/glimpse",
    # VLA CMZ compact sources (Lu+ 2019, ApJS 244, 35, table 3)
    "VLA_Lu2019":  "J/ApJS/244/35/table3",
    # Chandra Source Catalog 2.1 (Evans+ 2024)
    "Chandra_CSC": "IX/70/csc21mas",
    # Chandra GC X-ray catalog (Muno+ 2009, ApJS 181, 110)
    "Muno_Xray":   "J/ApJS/181/110/catalog",
}

# Columns we want back from each Vizier catalog.  Keeping the list short
# reduces bandwidth and makes the output table manageable.
VIZIER_COLUMNS = {
    "2MASS":       ["2MASS", "RAJ2000", "DEJ2000", "Jmag", "Hmag", "Kmag", "Qflg"],
    "GLIMPSE":     ["GLIMPSE", "RAJ2000", "DEJ2000",
                    "3.6mag", "4.5mag", "5.8mag", "8.0mag"],
    "VLA_Lu2019":  ["ID", "RAJ2000", "DEJ2000", "Ipeak", "Fint", "alpha"],
    "Chandra_CSC": ["2CXO", "RAICRS", "DEICRS"],
    "Muno_Xray":   ["CXOGC", "RAJ2000", "DEJ2000"],
}

# ---------------------------------------------------------------------------
# MeerKAT – no compact-source Vizier catalog exists for the Heywood+ 2022
# GC mosaic.  We therefore cross-match via SIMBAD radio associations, and
# the user can supplement this with a local MeerKAT catalog if available.
# We include the MALS continuum catalog (Wagenveld+ 2023) as a partial
# radio cross-check.
# ---------------------------------------------------------------------------
VIZIER_CATALOGS["MeerKAT"] = "J/A+A/673/A113/mals-ten"
VIZIER_COLUMNS["MeerKAT"] = ["*"]
MATCH_RADIUS["MeerKAT"] = 5.0


# ===================================================================== #
#                          Helper functions                              #
# ===================================================================== #

def load_aces_catalog() -> Table:
    """Load the ACES catalog and add an equatorial RA/Dec column pair."""
    cat = Table.read(CATALOG_PATH)
    print(f"Loaded {len(cat)} ACES sources from {CATALOG_PATH.name}")

    # Build SkyCoord in Galactic, convert to ICRS for crossmatching
    gc = SkyCoord(l=cat["GLON_peak"], b=cat["GLAT_peak"],
                  frame=Galactic, unit="deg")
    icrs = gc.icrs
    cat["RA_J2000"]  = icrs.ra.deg
    cat["Dec_J2000"] = icrs.dec.deg

    # Unique identifier for XMatch (must be integer-castable or string)
    if "index" not in cat.colnames:
        cat["index"] = np.arange(len(cat))
    return cat


def xmatch_vizier(aces_table: Table, vizier_id: str, radius_arcsec: float,
                  columns: list | None = None, label: str = "") -> Table:
    """
    Use the CDS XMatch service to cross-match the full ACES catalog
    against a Vizier catalog in one vectorised call.

    Parameters
    ----------
    aces_table : Table
        Must contain ``RA_J2000`` and ``Dec_J2000`` columns.
    vizier_id : str
        Vizier catalog identifier, e.g. ``"II/246/out"``.
    radius_arcsec : float
        Maximum match distance in arcseconds.
    columns : list or None
        Columns to request from the Vizier side (``None`` → all).
    label : str
        Human-readable label for progress messages.

    Returns
    -------
    result : `~astropy.table.Table`
        Matched rows, including the ``index`` column from the ACES table.
    """
    tag = label or vizier_id
    print(f"  XMatch: ACES ↔ {tag} (r = {radius_arcsec}″) …", end=" ", flush=True)

    # Build a minimal upload table
    upload = Table()
    upload["index"]    = aces_table["index"]
    upload["RA_J2000"] = aces_table["RA_J2000"]
    upload["Dec_J2000"] = aces_table["Dec_J2000"]

    # CDS XMatch wants a vizier:<catalog> string for cat2
    cat2 = f"vizier:{vizier_id}"

    try:
        result = XMatch.query(
            cat1=upload,
            cat2=cat2,
            max_distance=radius_arcsec * u.arcsec,
            colRA1="RA_J2000",
            colDec1="Dec_J2000",
        )
        print(f"{len(result)} matches")
    except Exception as exc:
        print(f"FAILED ({exc})")
        result = Table()

    return result


def query_simbad_bulk(aces_table: Table,
                      radius_arcsec: float = 3.0) -> Table:
    """
    Query SIMBAD for all ACES sources in one TAP/ADQL call.

    We upload a table of coordinates and use SIMBAD's ``query_tap`` with
    an ``UPLOAD`` join to do a cone-search per source on the server side.
    Falls back to a region query if TAP upload is not supported.

    Returns a table with columns:
        ``index``, ``simbad_name``, ``simbad_otype``, ``simbad_sep_arcsec``
    """
    print(f"  SIMBAD bulk query (r = {radius_arcsec}″) …", flush=True)

    simbad = Simbad()
    simbad.ROW_LIMIT = 1000 # chatbots like to set this to zero, which is _not_ unlimited, it literally means zero and will return nothing

    # ---- strategy: upload the whole catalog and run an ADQL cross-match
    # Note: column name "dec" is reserved in some ADQL parsers, so we use
    # "decl" for the uploaded declination column.
    upload = Table()
    upload["aces_idx"]  = aces_table["index"]
    upload["ra"]        = aces_table["RA_J2000"]
    upload["decl"]      = aces_table["Dec_J2000"]

    adql = f"""
    SELECT
        u.aces_idx,
        b.main_id   AS simbad_name,
        b.otype     AS simbad_otype,
        b.otype_txt AS simbad_otype_long,
        DISTANCE(
            POINT('ICRS', u.ra, u.decl),
            POINT('ICRS', b.ra, b.dec)
        ) AS sep_deg
    FROM TAP_UPLOAD.aces AS u
    JOIN basic AS b
      ON 1 = CONTAINS(
            POINT('ICRS', b.ra, b.dec),
            CIRCLE('ICRS', u.ra, u.decl, {radius_arcsec / 3600.0})
         )
    """

    try:
        result = simbad.query_tap(adql, aces=upload)
        if result is None or len(result) == 0:
            raise ValueError("empty SIMBAD TAP result")
        print(f"    → {len(result)} raw SIMBAD matches")
    except Exception as exc:
        print(f"    TAP upload failed ({exc}); falling back to XMatch/simbad …")
        # Fallback: CDS XMatch against "simbad"
        upload2 = Table()
        upload2["index"]    = aces_table["index"]
        upload2["RA_J2000"] = aces_table["RA_J2000"]
        upload2["Dec_J2000"] = aces_table["Dec_J2000"]
        try:
            result = XMatch.query(
                cat1=upload2,
                cat2="simbad",
                max_distance=radius_arcsec * u.arcsec,
                colRA1="RA_J2000",
                colDec1="Dec_J2000",
            )
            print(f"    → {len(result)} XMatch/simbad matches")
            # Rename columns to uniform names
            if "main_id" in result.colnames:
                result.rename_column("main_id", "simbad_name")
            if "main_type" in result.colnames:
                result.rename_column("main_type", "simbad_otype")
            if "angDist" in result.colnames:
                result["simbad_sep_arcsec"] = result["angDist"]
            result.rename_column("index", "aces_idx")
        except Exception as exc2:
            print(f"    XMatch fallback also failed ({exc2})")
            result = Table({"aces_idx": [], "simbad_name": [],
                            "simbad_otype": [], "simbad_sep_arcsec": []})
        return _deduplicate_simbad(result)

    # Convert separation from degrees to arcsec
    result["simbad_sep_arcsec"] = result["sep_deg"] * 3600.0
    result.remove_column("sep_deg")

    return _deduplicate_simbad(result)


def _deduplicate_simbad(tbl: Table) -> Table:
    """Keep only the closest SIMBAD match per ACES source."""
    if len(tbl) == 0:
        return tbl

    idx_col = "aces_idx"
    sep_col = "simbad_sep_arcsec"
    if idx_col not in tbl.colnames or sep_col not in tbl.colnames:
        return tbl

    # Group by ACES index, keep the row with smallest separation
    tbl.sort([idx_col, sep_col])
    _, unique_idx = np.unique(tbl[idx_col], return_index=True)
    return tbl[unique_idx]


# ===================================================================== #
#                      Classification logic                              #
# ===================================================================== #

def classify_sources(master: Table) -> Table:
    """
    Add an ``object_type`` column to the master table based on the
    multi-wavelength association columns.

    Priority order (highest to lowest):
      1. dust              (1 mm detection + α_1mm-3mm ≥ 2)
      2. radio+nir         (radio + bright NIR)
      3. evolved_star      (bright NIR, no radio)
      4. hii_region        (radio, no bright NIR)
      5. xray_source       (Chandra / Muno match)
      6. simbad:<otype>    (SIMBAD classification)
      7. unclassified      (no associations)
    """
    n = len(master)
    otypes = np.full(n, "unclassified", dtype="U40")

    # ----- helpers to test column presence and non-emptiness ----------
    def has_match(colname):
        """Boolean array: True where the column has a non-empty value."""
        if colname not in master.colnames:
            return np.zeros(n, dtype=bool)
        col = master[colname]
        if col.dtype.kind in ('U', 'S', 'O'):  # string-like
            return np.array([bool(v) and str(v).strip() not in ('', '--', 'N/A')
                             for v in col])
        else:
            return ~np.isnan(np.array(col, dtype=float))

    def finite_col(colname, default=np.nan):
        if colname not in master.colnames:
            return np.full(n, default)
        arr = np.array(master[colname], dtype=float)
        return arr

    # ----- flags -------------------------------------------------------
    has_cmzoom   = has_match("cmzoom_id") | (finite_col("cmzoom_flux") > 0)
    alpha_1mm    = finite_col("alpha_cmzoom_aces")
    is_dust      = has_cmzoom & np.isfinite(alpha_1mm) & (alpha_1mm >= 2.0)

    has_2mass    = has_match("match_2MASS")
    has_glimpse  = has_match("match_GLIMPSE")
    kmag         = finite_col("Kmag_2MASS")
    mag45        = finite_col("mag45_GLIMPSE")
    bright_nir   = ((has_2mass & np.isfinite(kmag) & (kmag < 10)) |
                    (has_glimpse & np.isfinite(mag45) & (mag45 < 9)))

    has_vla      = has_match("match_VLA_Lu2019")
    has_meerkat  = has_match("match_MeerKAT")
    has_radio    = has_vla | has_meerkat

    has_chandra  = has_match("match_Chandra_CSC")
    has_muno     = has_match("match_Muno_Xray")
    has_xray     = has_chandra | has_muno

    has_simbad   = has_match("simbad_name")

    # ----- assign types (lowest priority first, overwrite upward) ------
    # 6. SIMBAD fallback
    if "simbad_otype" in master.colnames:
        for i in range(n):
            if has_simbad[i]:
                ot = str(master["simbad_otype"][i]).strip()
                if ot and ot not in ('', '--'):
                    otypes[i] = f"simbad:{ot}"

    # 5. X-ray
    otypes[has_xray] = "xray_source"

    # 4. H II region (radio, no bright NIR)
    otypes[has_radio & ~bright_nir] = "hii_region"

    # 3. evolved star (bright NIR, no radio)
    otypes[bright_nir & ~has_radio] = "evolved_star"

    # 2. radio + NIR
    otypes[has_radio & bright_nir] = "radio+nir"

    # 1. dust-dominated
    otypes[is_dust] = "dust"

    master["object_type"] = otypes
    return master


# ===================================================================== #
#                             Main routine                               #
# ===================================================================== #

def main():
    # ------ load ACES catalog -----------------------------------------
    aces = load_aces_catalog()

    # ------ Vizier XMatch for each catalog ----------------------------
    print("\n=== Vizier crossmatches ===")
    xmatch_results = {}
    for cat_label, cat_id in VIZIER_CATALOGS.items():
        radius = MATCH_RADIUS.get(cat_label, 3.0)
        columns = VIZIER_COLUMNS.get(cat_label)
        xm = xmatch_vizier(aces, cat_id, radius, columns=columns,
                           label=cat_label)
        xmatch_results[cat_label] = xm

    # ------ SIMBAD bulk query -----------------------------------------
    print("\n=== SIMBAD query ===")
    simbad_matches = query_simbad_bulk(aces, radius_arcsec=3.0)

    # ------ Build master output table ---------------------------------
    print("\n=== Building master table ===")
    master = aces.copy()

    # --- merge SIMBAD -------------------------------------------------
    if len(simbad_matches) > 0:
        # Ensure consistent column name
        if "aces_idx" in simbad_matches.colnames:
            simbad_matches.rename_column("aces_idx", "index")
        keep_cols = [c for c in ("index", "simbad_name", "simbad_otype",
                                 "simbad_otype_long", "simbad_sep_arcsec")
                     if c in simbad_matches.colnames]
        simbad_slim = simbad_matches[keep_cols]
        master = join(master, simbad_slim, keys="index", join_type="left")
        print(f"  SIMBAD: {np.sum(~master['simbad_name'].mask) if hasattr(master['simbad_name'], 'mask') else 'all'} matches merged")
    else:
        master["simbad_name"]  = ""
        master["simbad_otype"] = ""

    # --- merge each Vizier catalog ------------------------------------
    for cat_label, xm in xmatch_results.items():
        match_col = f"match_{cat_label}"
        sep_col   = f"sep_{cat_label}"

        if len(xm) == 0:
            master[match_col] = ""
            master[sep_col]   = np.nan
            print(f"  {cat_label}: 0 matches")
            continue

        # Identify the "name/ID" column from the Vizier side
        name_col = _pick_name_column(xm, cat_label)

        # De-duplicate: keep closest match per ACES source
        if "angDist" in xm.colnames:
            xm.sort(["index", "angDist"])
        _, uidx = np.unique(xm["index"], return_index=True)
        xm = xm[uidx]

        # Build a slim table for joining
        slim = Table()
        slim["index"] = xm["index"]
        slim[match_col] = xm[name_col] if name_col else ""
        if "angDist" in xm.colnames:
            slim[sep_col] = xm["angDist"]
        else:
            slim[sep_col] = np.nan

        # Carry across photometric columns of interest
        extra_cols = _extra_photometry_cols(xm, cat_label)
        for src_col, dst_col in extra_cols:
            if src_col in xm.colnames:
                slim[dst_col] = xm[src_col]

        master = join(master, slim, keys="index", join_type="left")
        n_matched = np.sum(slim[match_col] != "") if slim[match_col].dtype.kind in ('U','S','O') else len(slim)
        print(f"  {cat_label}: {n_matched} matches merged")

    # --- classify -----------------------------------------------------
    print("\n=== Classifying sources ===")
    master = classify_sources(master)

    # Print summary
    types, counts = np.unique(master["object_type"], return_counts=True)
    print("\nClassification summary:")
    for t, c in sorted(zip(types, counts), key=lambda x: -x[1]):
        print(f"  {t:30s}  {c:5d}")

    # --- write output --------------------------------------------------
    out_fits = OUTPUT_DIR / "aces_crossmatched_catalog.fits"
    out_ecsv = OUTPUT_DIR / "aces_crossmatched_catalog.ecsv"

    # Clean masked columns for FITS compatibility
    for col in master.colnames:
        if hasattr(master[col], "filled"):
            if master[col].dtype.kind in ("U", "S", "O"):
                master[col] = master[col].filled("")
            elif master[col].dtype.kind in ("i", "u"):  # integer types
                master[col] = master[col].filled(-1)
            else:
                master[col] = master[col].filled(np.nan)
        # Convert object-dtype columns (e.g. from SIMBAD TAP) to strings
        if master[col].dtype.kind == "O":
            master[col] = Column([str(v) if v is not None else ""
                                  for v in master[col]], name=col)

    master.write(out_fits, overwrite=True)
    master.write(out_ecsv, format="ascii.ecsv", overwrite=True)
    print(f"\nWrote {len(master)} rows to:")
    print(f"  {out_fits}")
    print(f"  {out_ecsv}")


# ----- small helpers for column name mapping ---------------------------

def _pick_name_column(xm: Table, label: str) -> str | None:
    """Return the Vizier-side name/ID column for a given catalog."""
    preference = {
        "2MASS":       "2MASS",
        "GLIMPSE":     "GLIMPSE",
        "VLA_Lu2019":  "ID",
        "Chandra_CSC": "2CXO",
        "Muno_Xray":   "CXOGC",
        "MeerKAT":     "MALS",
    }
    col = preference.get(label)
    if col and col in xm.colnames:
        return col
    # Fallback: first string-like column that is not RA/Dec/index
    for c in xm.colnames:
        if c.lower() in ("index", "ra_j2000", "dec_j2000", "angdist",
                          "raj2000", "dej2000", "raicrs", "deicrs"):
            continue
        if xm[c].dtype.kind in ("U", "S", "O"):
            return c
    return None


def _extra_photometry_cols(xm: Table, label: str) -> list[tuple[str, str]]:
    """Return (src_col, dest_col) pairs for photometry we want to keep."""
    mapping = {
        "2MASS":   [("Jmag", "Jmag_2MASS"), ("Hmag", "Hmag_2MASS"),
                    ("Kmag", "Kmag_2MASS"), ("Qflg", "Qflg_2MASS")],
        "GLIMPSE": [("3.6mag", "mag36_GLIMPSE"), ("4.5mag", "mag45_GLIMPSE"),
                    ("5.8mag", "mag58_GLIMPSE"), ("8.0mag", "mag80_GLIMPSE")],
        "VLA_Lu2019": [("Ipeak", "Ipeak_VLA"), ("Fint", "Fint_VLA"),
                       ("alpha", "alpha_VLA")],
    }
    return mapping.get(label, [])


if __name__ == "__main__":
    main()
