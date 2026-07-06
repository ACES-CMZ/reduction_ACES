# Giant-cube moments: adding a line & mask-guided moments

How the ACES giant-cube moment pipeline works, how to add a new line, and how to
make moment maps of a faint line using a brighter molecule's signal mask as a
guide.

Relevant code:
- `aces/analysis/giantcube_cuts.py` — moment / statistics computation
  (`aces_giantcube_analysis` entry point, `main()` driven by env var `MOLNAME`).
- `aces/imaging/mosaic_12m.py` — per-line cube builders (`make_giant_mosaic_cube_*`).
- `aces/hipergator_scripts/slurm_*_mosaic_jobs.sh` — build the mosaics.
- `aces/hipergator_scripts/run_allgiantcube_analysis.sh` — run the moments.

Data locations (`basepath = /orange/adamginsburg/ACES/`):
- Cubes: `{basepath}/mosaics/cubes/{MOLNAME}_CubeMosaic.fits`
  (+ `_downsampled9`, `_spectrally` companions).
- Moment products & masks: `{basepath}/mosaics/cubes/moments/`.
- Feathered per-field input cubes: `{basepath}/upload/Feather_12m_7m_TP/SPW{NN}/cubes/`.

---

## 1. Adding a new line to take moments of

Three edits + a two-stage mosaic job. Worked example: **CH3CHO 3(-1,3)-2(0,2) E at
101.343448 GHz** (Cont2 = SPW35), added as `CH3CHO_3m13` (kept distinct from the
existing SPW33 `CH3CHO` cube at 98.900951 GHz).

### 1a. Cube builder — `aces/imaging/mosaic_12m.py`
Add `make_giant_mosaic_cube_<name>(**kwargs)`, cloned from a sibling in the same
SPW (e.g. `make_giant_mosaic_cube_nsplus` for SPW35). Set:
- `filelist` glob for the correct `SPW{NN}` feather directory.
- `weightfilelist = [get_weightfile(fn, spw=NN) for fn in filelist]`.
- `restfrq` (Hz), `cdelt_kms` (channel width), `cubename`, `nchan`,
  `beam_threshold`, `channelmosaic_directory`.
- Keep the trailing `make_downsampled_cube(...)` (and `downsample_spectrally(...)`
  if the cube is large).

The `cubename` **must be an exact, unique string** — the velocity-mask
special-casing in `do_all_stats` (giantcube_cuts.py) matches molecule names by
exact tuple membership, so a near-duplicate name will not accidentally inherit
another line's velocity window.

Band-edge note: check that `nchan` × `cdelt_kms` centred on `restfrq` actually
covers the CMZ velocity range (roughly -200 to +200 km/s) without running off the
SPW edge. 101.343 GHz sits ~0.09 GHz from the SPW35 upper edge (≈ -280 km/s), so
the blue wing is bounded; `nchan=350` at 1.47 km/s (~515 km/s span) is adequate.

### 1b. Mosaic launcher — `aces/hipergator_scripts/slurm_<name>_mosaic_jobs.sh`
Clone a sibling (e.g. `slurm_nsplus_mosaic_jobs.sh`). Two stages:
1. `--array=0-49` job calling `make_giant_mosaic_cube_<name>(channels='slurm', skip_final_combination=True)`.
2. `--dependency=afterok:$jobid` merge job calling `make_giant_mosaic_cube_<name>(channels='all', skip_channel_mosaicing=True)`.

Submit with `bash slurm_<name>_mosaic_jobs.sh`. Produces `{cubename}_CubeMosaic.fits`.

### 1c. Register the moments run — `run_allgiantcube_analysis.sh`
Add `<cubename>` to the `for MOLNAME in ...` loop. The moment runner picks the
cube up by name; it gates on `{MOLNAME}_CubeMosaic.fits` existing, so run 1b first.
No change to `giantcube_cuts.py` is required for a plain line.

### 1d. (Optional) velocity window
If the line needs a non-default velocity window, add its exact `cubename` to the
appropriate branch of the `if molname in (...)` block in `do_all_stats`
(giantcube_cuts.py). Otherwise it uses `velocity_mask(cube)`.

---

## 2. Mask-guided moments (brighter molecule as guide)

Make moment maps of a **faint** line using a **brighter** molecule's
velocity-resolved (3D PPV) signal mask as the guide. This runs as the **last**
step of `main()` so it reuses the guide molecule's already-persisted mask.

### How it works
Each molecule's giant-cube run persists its full 3D signal mask to
`{moments}/{MOLNAME}_CubeMosaic_signal_mask_pruned.fits`. `get_guide_mask()`:
1. loads the guide molecule's mask cube,
2. asserts the shared celestial grid (ACES line mosaics share spatial WCS/shape —
   no spatial reprojection),
3. spectrally interpolates it onto the target cube's velocity axis
   (`SpectralCube.spectral_interpolate`), matching by **LSRK velocity** (both cubes
   are on km/s axes), and
4. returns a `BooleanArrayMask` (threshold > 0.5).

`do_guided_moments()` applies `velocity_mask` then the guide mask and writes:
- `{MOLNAME}_CubeMosaic_guidedby_{GUIDE}_mom0.fits` / `.png`
- `{MOLNAME}_CubeMosaic_guidedby_{GUIDE}_mom1.fits` / `.png`

### How to run
Set the `GUIDE_MOLNAME` env var alongside `MOLNAME`:

```bash
MOLNAME=CH3CHO_3m13 GUIDE_MOLNAME=CS21 aces_giantcube_analysis
```

Leave `GUIDE_MOLNAME` unset (or `False`) for a normal self-masked run — default
behavior is unchanged. In `run_allgiantcube_analysis.sh`, per-molecule guides are
set in the `case "$MOLNAME"` block and threaded through `sbatch --export`.

Requirements:
- The guide molecule must already have its `*_signal_mask_pruned.fits` on disk
  (run its giant-cube analysis first). Guides available today include CS21, HNCO_7m12mTP,
  HCOP, SO21, SO32, H13CN, H13COp, NSplus, SiO21, HC3N, HC15N, HN13C, H40a, CH3CHO.
- Guide and target must share the celestial grid (all ACES line mosaics do).
