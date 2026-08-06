# Experimental zarr giant-cube combiner

Single-pass alternative to the channel-by-channel `make_giant_mosaic_cube`.
Hands the full list of field cubes to one
`reproject.mosaicking.reproject_and_coadd(return_type='zarr')` call, which
computes the 3D mosaic block-by-block into a zarr store (bounded memory, truly
parallel I/O) and returns dask arrays that are then streamed to FITS. Based on
astrofrog's recipe (gist `81729066204b0332735b9461d015ee7c`) and
astropy/reproject#621.

## Isolation

Everything runs in a dedicated conda env and editable worktrees so the
production `python312` env and its running jobs are never disturbed:

| component      | worktree                                                      | branch                     |
|----------------|---------------------------------------------------------------|----------------------------|
| reproject      | `/blue/adamginsburg/adamginsburg/repos/reproject-zarr`        | `coadd-zarr-return-wt`     |
| spectral-cube  | `/blue/adamginsburg/adamginsburg/repos/spectral-cube-zarr`    | `zarr-combine`             |
| reduction_ACES | `/blue/adamginsburg/adamginsburg/repos/reduction_ACES-zarr`   | `giantcube-zarr-combiner`  |
| conda env      | `python312-zarr` (clone of `python312`, editables repointed)  |                            |

## Run

```bash
MOLNAME=CH3CHO_3m13 sbatch aces/hipergator_scripts/slurm_zarr_combine.sh
```

`MOLCONFIG` in `run_zarr_combine.py` carries per-molecule params (spw, restfreq,
cdelt, nchan, beam threshold) copied from the `make_giant_mosaic_cube_*`
builders. Add a molecule by adding a row.

Key settings (overridable via env): `PARALLEL=16`, `block_size=(1, 2048, -1)`,
`zarr_batch_size=64`.

## Open items to tune during the first real test

- **Beam matching**: the channel path convolves each field to a common beam
  before coadd; this combiner does *not* yet (the 3D convolve is expensive).
  Decide whether per-field beam differences matter for the target line.
- **Array extraction**: cubes are passed as `(cube.unitless_filled_data[:],
  cube.wcs)`; weights as `nan_to_num(wc.unitless_filled_data[:])`. Confirm the
  dask arrays stay lazy end-to-end (no accidental `.compute()`).
- **block_size last dim** should match the output width for best efficiency
  (`-1` = full width); raise `zarr_batch_size` on many cores.
- `zarr_path` must not pre-exist (reproject enforces); `OVERWRITE=True` removes it.
