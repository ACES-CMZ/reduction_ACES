"""
Prepare (normalize + bundle) the maser catalogs used by crossmatch_masers.py.

Each source catalog is fetched (from VizieR or a local file), reduced to a
uniform schema, and written to ``aces/data/masers/<label>.fits``:

    ra_deg, dec_deg, vlsr [km/s], flux, flux_unit, name, species, catalog, reference

Bundled catalogs (covering the ACES / CMZ footprint):

* H2O   -- Walsh et al. 2014 (SWAG/ATCA, all-plane; local file)
* H2O   -- Lu et al. 2019 CMZ water masers (VizieR J/ApJS/244/35 table4)
* CH3OH -- Cotton & Yusef-Zadeh 2016 (VLA 36 GHz class I, CMZ; VizieR J/ApJS/227/10)
          -- the dense CMZ methanol catalog (2240 maser spots)
* CH3OH -- GLOSTAR 6.7 GHz class II (all-plane; local file)

Re-run this (needs network for the VizieR catalogs) to regenerate the bundled
FITS files.  Run: ``python -m aces.analysis.prepare_maser_catalogs``.
"""
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
from astropy.table import Table

MASER_DATA_DIR = Path(__file__).parent.parent / "data" / "masers"

# Local raw catalogs (not in the repo)
WALSH2014 = "/orange/adamginsburg/galactic_plane_surveys/water_masers/walsh2014_atca_water_masers.fits"
GLOSTAR = "/orange/adamginsburg/salt/survey_2026/data/glostar_methanol.fits"


def _sexagesimal(ra, dec):
    return SkyCoord(np.array([str(x).strip() for x in ra]),
                    np.array([str(x).strip() for x in dec]),
                    unit=(u.hourangle, u.deg))


def _galactic(lon, lat):
    return SkyCoord(l=np.asarray(lon, float) * u.deg, b=np.asarray(lat, float) * u.deg,
                    frame=Galactic).icrs


def _normalized(coords, vlsr, flux, flux_unit, name, species, catalog, reference,
                maser_class=""):
    n = len(coords)
    # maser_type distinguishes CH3OH class I vs II (critical: different pumping /
    # physics); H2O has no class.
    if species == "CH3OH" and maser_class:
        maser_type = f"CH3OH_class{maser_class}"
    else:
        maser_type = species
    out = Table()
    out["ra_deg"] = coords.ra.deg
    out["dec_deg"] = coords.dec.deg
    out["vlsr"] = np.asarray(vlsr, float) if vlsr is not None else np.full(n, np.nan)
    out["flux"] = np.asarray(flux, float) if flux is not None else np.full(n, np.nan)
    out["flux_unit"] = flux_unit
    out["name"] = np.array([str(x).strip() for x in name]) if name is not None else np.array([""] * n)
    out["species"] = species
    out["maser_class"] = maser_class
    out["maser_type"] = maser_type
    out["catalog"] = catalog
    out["reference"] = reference
    return out


def prepare():
    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    MASER_DATA_DIR.mkdir(parents=True, exist_ok=True)

    # --- H2O: Walsh 2014 (local, sexagesimal) ---
    w = Table.read(WALSH2014)
    good = np.array([str(r).strip() not in ("", "--", "nan") for r in w["RAJ2000"]])
    w = w[good]
    t = _normalized(_sexagesimal(w["RAJ2000"], w["DEJ2000"]), w["Vp"], w["Sp"], "Jy",
                    w["Name"], "H2O", "H2O_Walsh2014",
                    "Walsh et al. 2014, MNRAS 442, 2240 (SWAG/ATCA water masers)")
    t.write(MASER_DATA_DIR / "h2o_walsh2014.fits", overwrite=True)
    print(f"h2o_walsh2014.fits: {len(t)}")
    # NOTE: the SWAG water masers (Ward et al., in prep.) are NOT bundled; they are
    # read from disk at runtime by crossmatch_masers._load_swag_masers().

    # --- H2O: Lu 2019 CMZ water masers (VizieR J/ApJS/244/35 table4, sexagesimal) ---
    lu = Vizier.get_catalogs("J/ApJS/244/35/table4")[0]
    flux = lu["Fpk"] if "Fpk" in lu.colnames else None
    t = _normalized(_sexagesimal(lu["RAJ2000"], lu["DEJ2000"]), None, flux, "Jy",
                    lu["ID"], "H2O", "H2O_Lu2019_CMZ",
                    "Lu et al. 2019, ApJS 244, 35 table4 (CMZ water masers)")
    t.write(MASER_DATA_DIR / "h2o_lu2019_cmz.fits", overwrite=True)
    print(f"h2o_lu2019_cmz.fits: {len(t)}")

    # --- CH3OH: Cotton & Yusef-Zadeh 2016 (VizieR J/ApJS/227/10, Galactic; 36 GHz class I) ---
    co = Vizier.get_catalogs("J/ApJS/227/10")[0]
    flux_jy = np.asarray(co["Flux"], float) / 1000.0  # mJy -> Jy
    t = _normalized(_galactic(co["GLON"], co["GLAT"]), co["Vel"], flux_jy, "Jy",
                    co["Name"], "CH3OH", "CH3OH_Cotton2016_CMZ",
                    "Cotton & Yusef-Zadeh 2016, ApJS 227, 10 (VLA 36 GHz class I CH3OH, CMZ)",
                    maser_class="I")
    t.write(MASER_DATA_DIR / "ch3oh_cotton2016_cmz.fits", overwrite=True)
    print(f"ch3oh_cotton2016_cmz.fits: {len(t)}")

    # --- CH3OH: GLOSTAR 6.7 GHz class II (local, sexagesimal) ---
    g = Table.read(GLOSTAR)
    good = np.array([str(r).strip() not in ("", "--", "nan") for r in g["RAJ2000"]])
    g = g[good]
    t = _normalized(_sexagesimal(g["RAJ2000"], g["DEJ2000"]), g["Vlsr"], g["Svp"], "Jy/beam",
                    g["Name"], "CH3OH", "CH3OH_GLOSTAR",
                    "Nguyen et al. 2022, A&A 666, A59 (GLOSTAR 6.7 GHz class II CH3OH)",
                    maser_class="II")
    t.write(MASER_DATA_DIR / "ch3oh_glostar.fits", overwrite=True)
    print(f"ch3oh_glostar.fits: {len(t)}")


if __name__ == "__main__":
    prepare()
