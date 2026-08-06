#!/bin/bash
# Experimental single-pass giant-cube combiner (reproject return_type='zarr').
#
# Runs in the dedicated python312-zarr env (editable worktrees of
# reproject@coadd-zarr-return, spectral-cube@zarr-combine, reduction_ACES@
# giantcube-zarr-combiner) so the production python312 env / running jobs are
# untouched.
#
# Usage:
#   MOLNAME=CH3CHO_3m13 sbatch aces/hipergator_scripts/slurm_zarr_combine.sh
# Optional env: PARALLEL, ZARR_BATCH_SIZE, BLOCK_SPATIAL, ZARR_PATH, OUTPUT_FILE, OVERWRITE
#
# ntasks should match PARALLEL (reproject uses that many threads).

#SBATCH --job-name=aces_zarr_combine
#SBATCH --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_zarr_combine_%j.log
#SBATCH --account=astronomy-dept
#SBATCH --qos=astronomy-dept-b
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=256gb
#SBATCH --time=96:00:00

export MOLNAME=${MOLNAME:?set MOLNAME (e.g. CH3CHO_3m13)}
export PARALLEL=${PARALLEL:-16}

PY=/blue/adamginsburg/adamginsburg/miniconda3/envs/python312-zarr/bin/python
echo "MOLNAME=$MOLNAME PARALLEL=$PARALLEL"
$PY -c "import reproject, spectral_cube, aces; print('reproject', reproject.__version__); print('spectral_cube', spectral_cube.__version__)"
$PY -m aces.imaging.run_zarr_combine
