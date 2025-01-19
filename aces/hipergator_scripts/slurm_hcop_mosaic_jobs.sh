
jobid=$(sbatch --job-name=aces_hcop_mos_arr \
    --output=/red/adamginsburg/ACES/logs/aces_hcop_mosaic_%j_%A_%a.log  \
    --array=0-99 \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --ntasks=8 --nodes=1 --mem=64gb --time=96:00:00 --parsable \
    --wrap "/red/adamginsburg/miniconda3/envs/python310/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_hcop; make_giant_mosaic_cube_hcop(channels='slurm', skip_final_combination=True, verbose=True,)\"")

echo "Job IDs are ${jobid}"

sbatch --job-name=aces_hcop_mosaic_merge \
    --output=/red/adamginsburg/ACES/logs/aces_hcop_mosaic_merge_%j.log  \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --dependency=afterok:$jobid \
    --ntasks=16 --nodes=1 --mem=64gb --time=96:00:00 \
    --wrap "/red/adamginsburg/miniconda3/envs/python310/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_hcop; make_giant_mosaic_cube_hcop(channels='all', skip_channel_mosaicing=True, verbose=True,)\""
