
jobid=$(sbatch --job-name=aces_hcop_mopra_mos_arr \
    --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_hcop_mopra_mosaic_%j_%A_%a.log  \
    --array=0-54 \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --ntasks=8 --nodes=1 --mem=64gb --time=96:00:00 --parsable \
    --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python310/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_hcop_mopra; make_giant_mosaic_cube_hcop_mopra(channels='slurm', skip_final_combination=True, verbose=True,)\"")

echo "Job IDs are ${jobid}"

sbatch --job-name=aces_hcop_mopra_mosaic_merge \
    --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_hcop_mopra_mosaic_merge_%j.log  \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --dependency=afterok:$jobid \
    --ntasks=16 --nodes=1 --mem=64gb --time=96:00:00 \
    --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python310/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_hcop_mopra; make_giant_mosaic_cube_hcop_mopra(channels='all', skip_channel_mosaicing=True, verbose=True,)\""
