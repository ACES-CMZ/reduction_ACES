#!/bin/bash
# Full-band, 5" pixel, 10" beam per-spw overview cubes via the zarr combiner.
# One array task per spw (25, 27, 33, 35).
#
# Runs in the python312-zarr env (editable worktrees of reproject@zarr,
# spectral-cube, reduction_ACES) so the reproject return_type='zarr' feature is
# available and the production env is untouched.  See README_zarr_combine.md.
#
# Usage:
#   sbatch aces/hipergator_scripts/slurm_fullband_downsampled.sh
# (or a single spw:  SPW=27 sbatch --array=0 ... )

#SBATCH --job-name=aces_fullband_ds
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_fullband_ds_%A_%a.log
#SBATCH --account=astronomy-dept
#SBATCH --qos=astronomy-dept-b
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --time=48:00:00
#SBATCH --array=0-3

SPWS=(25 27 33 35)
export SPW=${SPW:-${SPWS[$SLURM_ARRAY_TASK_ID]}}
export PARALLEL=${PARALLEL:-16}

PY=/blue/adamginsburg/adamginsburg/miniconda3/envs/python312-zarr/bin/python
echo "SPW=$SPW PARALLEL=$PARALLEL"
$PY -c "import reproject, spectral_cube; print('reproject', reproject.__version__)"
$PY -m aces.imaging.run_fullband_downsampled
