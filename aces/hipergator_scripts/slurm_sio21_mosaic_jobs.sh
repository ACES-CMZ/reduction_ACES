

jobid=$(sbatch --job-name=aces_sio21_mos_arr \
    --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_sio21_mosaic_%a_%A_%j.log  \
    --array=0-99 \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --ntasks=8 --nodes=1 --mem=32gb --time=96:00:00 --parsable \
    --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python39/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_sio21; make_giant_mosaic_cube_sio21(channels='slurm', skip_final_combination=True, verbose=True,)\"")

echo "Job IDs are ${jobid}"

sbatch --job-name=aces_sio21_mosaic_merge \
    --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_sio21_mosaic_merge_%j.log  \
    --dependency=afterok:$jobid \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --ntasks=8 --nodes=1 --mem=32gb --time=96:00:00 \
    --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python39/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_sio21; make_giant_mosaic_cube_sio21(channels='all', skip_channel_mosaicing=True, verbose=True,)\""
