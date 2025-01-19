jobid=$(sbatch --job-name=aces_h40a_mos_arr \
    --output=/red/adamginsburg/ACES/logs/aces_h40a_mosaic_%j_%A_%a.log  \
    --array=0-49 \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --ntasks=8 --nodes=1 --mem=64gb --time=96:00:00 --parsable \
    --wrap "/red/adamginsburg/miniconda3/envs/python310/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_h40a; make_giant_mosaic_cube_h40a(channels='slurm', skip_final_combination=True, verbose=True,)\"")

echo "Job IDs are ${jobid}"

sbatch --job-name=aces_h40a_mosaic_merge \
    --output=/red/adamginsburg/ACES/logs/aces_h40a_mosaic_merge_%j.log  \
    --dependency=afterok:$jobid \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --ntasks=8 --nodes=1 --mem=32gb --time=96:00:00 \
    --wrap "/red/adamginsburg/miniconda3/envs/python310/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_h40a; make_giant_mosaic_cube_h40a(channels='all', skip_channel_mosaicing=True, verbose=True,)\""
