
jobid=$(sbatch --job-name=aces_h13cn_mos_arr \
    --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_h13cn_mosaic_%j_%A_%a.log  \
    --array=0-49 \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --ntasks=8 --nodes=1 --mem=64gb --time=96:00:00 --parsable \
    --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python312/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_h13cn; make_giant_mosaic_cube_h13cn(channels='slurm', skip_final_combination=True, verbose=True,)\"")

echo "Job IDs are ${jobid}"

sbatch --job-name=aces_h13cn_mosaic_merge \
    --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_h13cn_mosaic_merge_%j.log  \
    --dependency=afterok:$jobid \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --ntasks=8 --nodes=1 --mem=32gb --time=96:00:00 \
    --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python312/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_h13cn; make_giant_mosaic_cube_h13cn(channels='all', skip_channel_mosaicing=True, verbose=True,)\""
