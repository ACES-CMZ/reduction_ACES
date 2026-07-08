# Line Detection and Catalog Augmentation

`line_detection_and_catalog_augmentation.py` searches for spectral-line emission
toward the compact continuum sources of the ACES catalog and records the fitted
line parameters as new catalog columns.

## Procedure

The input catalog is the ACES compact continuum source catalog
(`aces_compact_catalog_v0.fits`). The one-dimensional spectra it operates on are
produced beforehand by `spectral_extraction_everywhere.py`. For each catalog
source, that script defines an elliptical aperture from the catalog's fit
parameters — `fitted_major` and `fitted_minor` as the ellipse width and height
in arcseconds, `pa` as the position angle, centered on `GLON_peak`/`GLAT_peak` —
cuts the sub-cube covered by that ellipse from the source's 12m image cube, and
averages over the spatial axes to yield one spectrum per spectral window.  The
output files are tagged `ellipseaverage` and carry the source geometry in their
headers (`CATINDX`, `CATGLON`, `CATGLAT`, `CATMAJS`, `CATMINS`, `CATPA`). 

The script also writes a background spectrum per source, the mean over a
circular annulus centered on the source (inner radius 1.5× and outer radius 3×
the source's larger semi-axis), used later to test whether a detected line is
source-associated.  `line_detection_and_catalog_augmentation.py` consumes those
spectra for the four line-bearing 12m spectral windows — SPW 25, 27, 33, and 35;
SPW 29 (HCO+) and 31 (HNCO) are excluded. The line list, giving rest frequencies
and a star-based importance ranking, is read from `linelist.csv` combined with
`extended_linelist.csv`.

Fitting is done on the extracted source spectrum itself, not on
background-subtracted data; the annulus background is used only in the separate
verification step below. For each source and each line whose rest frequency
falls within the observed band, the script examines the spectrum in a ±200 km/s
window about the rest frequency. A constant baseline and a noise level are
estimated locally from that window by sigma-clipping (3σ, MAD-based), and it is
this per-window constant that is subtracted before the Gaussian fit — the fitted
amplitude is the line peak above that local baseline. The brightest channel is
located, preferring peaks within ±150 km/s of the rest velocity, and a single
Gaussian plus constant baseline is fit. A candidate is retained only if the peak
exceeds 4σ before fitting and the fitted amplitude exceeds 4.5σ of the fit
residuals; the fitted centroid must lie within the search window and the fitted
FWHM must fall between 0.5 and 300 km/s.  Where several lines crowd one window,
the fit is iterated — fit the strongest peak, subtract it, refit the residual —
up to five times, so that blended components are recovered separately.

Three consistency filters are then applied per source. First, when one line
yields several candidate peaks, the candidate whose velocity is closest to the
median velocity of the source's unambiguous single-peak detections is chosen.
For example, if a source's single-peak lines (CS, H13CN, HN13C) agree on
v ≈ +30 km/s, and SiO 2-1 returns two candidates at −80 and +25 km/s, the
+25 km/s candidate is kept and the −80 km/s one — likely a blend or noise peak
elsewhere in the window — is discarded.
Second, when two lines fit centroids within one FWHM of each other — i.e. the
same physical peak — only the higher-importance line (by the line list's star
ranking, breaking ties by signal-to-noise) is kept. Third, if a source has at
least three detected lines and one line's velocity is a strong outlier (beyond
5σ of the MAD-based scatter and more than 5 km/s from the median), that line is
discarded as a likely misidentification.

Surviving detections are verified against the local background. A mean spectrum
is extracted from an annulus around each source (inner radius 1.5× the source
semi-major axis, outer radius 3× that), taken from the parent cube. A detection
is confirmed only if the source-aperture line amplitude exceeds the background
peak in the same window by a factor of two, so that only source-associated
emission is kept.

Verified detections are tallied across the whole catalog. A line detected in at
least five sources is judged common and receives dedicated catalog columns;
lines detected in fewer sources are judged rare and are written to a JSON side
file. For common lines, the columns are populated for detected sources from the
Gaussian fit, and for non-detected sources by measuring the noise and the
in-window peak (an upper limit) directly from the spectrum, so that every source
carries a homogeneous noise characterization for every common line.

## Line importance ranking

The line list's `col9` field denotes priority with a string of asterisks:
`***` (highest), `**`, `*`, and blank (unranked, treated as lowest). This
ranking is used only to break ties when two lines fit the same peak — the
higher-tier line is kept. The tiers are:

- **`***` (4 lines):** HNCO 4-3 (87.925 GHz, SPW31), HCO+ 1-0 (89.189 GHz,
  SPW29), CS 2-1 (97.981 GHz, SPW33), HC3N 11-10 (100.076 GHz, SPW35). Two of
  these — HNCO and HCO+ — fall in SPW 31 and 29, which are not searched, so only
  CS 2-1 and HC3N 11-10 are effectively top-tier targets.
- **`**` (10 lines):** HC15N 1-0, SO 2(2)-1(1), SiO 2-1 v=1 maser, H13CN 1-0
  (SPW25); H13CO+ 1-0, SiO 2-1, HN13C 1-0 (SPW27); CH3CHO 5(1,4)-4(1,3) A, H40α,
  SO 3(2)-2(1) (SPW33).
- **`*` (8 lines):** H50β, HC3N 11-10 v7=1, NH2CN 5(1,4)-4(1,3), CH3SH 4(0)-3(0)
  A and E, CH3OH 7(2,5)-6(3,4) A, CH3CH2CN 10(1,10)-9(1,9), HC5N J=38-37.
- **blank (54 lines):** all remaining transitions in `linelist.csv` and
  `extended_linelist.csv`.

Every in-band line is fitted regardless of tier; the ranking affects only
deduplication, not whether a line is searched.

## Output files

Paths are relative to `conf.basepath` (`/orange/adamginsburg/ACES/`).

- `catalogs/aces_compact_catalog_v0_withLines.fits` and
  `catalogs/aces_compact_catalog_v0_withLines.ecsv` — the input catalog with
  the added line columns, written in both FITS and ECSV form.
- `catalogs/aces_compact_catalog_v0_rare_line_detections.json` — the rare-line
  detections that did not earn catalog columns.
- `spectra/plots/top40/source_NNNNN_spectra.png` — per-source spectral plots
  (four SPW panels, data, background, and fitted model) for the forty sources
  with the most distinct lines; `NNNNN` is the zero-padded source index.
- `spectra/plots/not_top40/source_NNNNN_spectra.png` — the same plots for all
  remaining sources, produced when `--all-sources` is given.
- `spectra/spatial_plots/spatial_<line>.png` — sky distribution of the sources
  detecting a given line.
- `spectra/spatial_plots/spatial_by_n_lines.png` — sky distribution with symbol
  size scaled by the number of lines detected per source.
- `spectra/spatial_plots/spatial_by_mean_velocity_nlines_gt4.png` — sky
  distribution of line-rich sources (more than four lines) colored by mean
  line-of-sight velocity.

Intermediate results are cached under `reduction_ACES/aces/analysis/.cache/`
(`all_detections.pkl`, `line_counts.pkl`, `verified_detections.pkl`, plus a
`simbad_queries.json` query cache). The cache is invalidated automatically when
the input catalog, the line list, or this script changes, tracked by MD5 hash in
`cache_metadata.json`.

## Catalog metadata

For each common line (identified by a filesystem-safe version of its name), six
columns are added, in units of Jy/beam for fluxes and km/s for velocities:

- `<line>_peak` — line amplitude. For detections this is the fitted Gaussian
  amplitude; for non-detections it is the maximum flux in the search window,
  serving as an upper limit.
- `<line>_rms` — noise level. For detections it is the RMS of the fit residuals;
  for non-detections it is the sigma-clipped RMS of the search window. Populated
  for every source, so a signal-to-noise ratio is recoverable as
  `<line>_peak / <line>_rms`.
- `<line>_fwhm_kms` — fitted line width; populated for detections only.
- `<line>_centroid_kms` — fitted centroid velocity; populated for detections only.
- `<line>_vmean` — center of the velocity window over which the flux and noise
  were measured: the line centroid for detections, 0 km/s for non-detections.
- `<line>_vwidth` — width of that measurement window: max(3× FWHM, 10 km/s) for
  detections, the full 200 km/s search range for non-detections.

Columns are left NaN where no measurement applies (for instance, width and
centroid at a non-detection).

## Rare-line JSON

`aces_compact_catalog_v0_rare_line_detections.json` is keyed by source index
(as a string), then by line name. Each entry records `peak`, `fwhm_kms`,
`centroid_kms`, `snr`, `line_name` (the original human-readable line name), and
`rest_freq_ghz`.

## Background-spectrum metadata

Background annulus spectra, extracted during verification, are cached next to the
source spectra with a `_background.fits` suffix. Their headers record the source
index (`CATINDX`), the source Galactic coordinates (`CATGLON`, `CATGLAT`), the
background type (`BKGTYPE = 'annulus'`), and the inner and outer radius scale
factors (`BKGINNER = 1.5`, `BKGOUTER = 3.0`).

## SIMBAD annotation

Plot titles and the printed top-40 ranking are annotated with the most-cited
SIMBAD object within 5″ of each source (name, object type, citation count).
Query results are cached to avoid repeated network calls.

## Invocation

Run with no arguments to detect, verify, augment, and plot using any valid
cache. `--regenerate-cache` forces a full rescan; `--skip-verification` omits the
background check; `--plot-only` regenerates plots from cached results; and
`--all-sources` plots every source with a detection rather than only the top 40.
