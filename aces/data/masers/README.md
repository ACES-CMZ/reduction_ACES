# Maser catalogs

External maser catalogs bundled for crossmatching against ACES continuum sources
(see `aces/analysis/crossmatch_masers.py`, entry point `aces_crossmatch_masers`).

| File | Species | Rows | Source |
|------|---------|------|--------|
| `walsh2014_atca_water_masers.fits` | H2O (water) | 2790 | Walsh et al. 2014, MNRAS 442, 2240 — ATCA water-maser survey of the Galactic plane / CMZ (the SWAG water-maser product). Columns: `Name`, `RAJ2000`, `DEJ2000` (sexagesimal), `Sp` [Jy] peak flux, `Vp` [km/s] peak velocity. |
| `glostar_methanol.fits` | CH3OH (6.7 GHz class II) | 2904 | GLOSTAR 6.7 GHz methanol maser catalogue, Nguyen et al. 2022, A&A 666, A59 (VizieR `J/A+A/666/A59`). Columns: `Name`, `RAJ2000`, `DEJ2000` (sexagesimal), `Vlsr` [km/s], `Svp` [Jy/beam] peak flux. |

The CMZ-specific class-I methanol maser survey of Lu et al. 2019 (ApJS 244, 35;
VizieR `J/ApJS/244/35`) is fetched from VizieR at runtime with the
`aces_crossmatch_masers --vizier` flag (not bundled).

These are published catalogs redistributed here for reproducibility; cite the
original references above when using them.
