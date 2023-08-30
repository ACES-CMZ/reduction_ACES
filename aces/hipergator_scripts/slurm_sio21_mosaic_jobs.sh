

incr=50
jobids=""

for ch in `seq 0 ${incr} 350`; do
    ch1=$ch
    ch2=$((${ch1} + ${incr}))
    jobid=$(sbatch --job-name=aces_sio21_mosaic_ch${ch1}to${ch2} \
        --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_sio21_mosaic_ch${ch1}to${ch2}_%j.log  \
        --account=astronomy-dept --qos=astronomy-dept-b \
        --ntasks=8 --nodes=1 --mem=32gb --time=96:00:00 \
        --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python39/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_sio21; make_giant_mosaic_cube_sio21(channels=range(${ch1},${ch2}, skip_final_combination=True, verbose=True,)\"")
    jobids+=$jobid+":"
done

echo "Job IDs are ${jobids}"

sbatch --job_name=aces_sio21_mosaic_merge \
    --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_sio21_mosaic_merge_%j.log  \
    --dependency=afterok:$jobids \
    --account=astronomy-dept --qos=astronomy-dept-b \
    --ntasks=8 --nodes=1 --mem=32gb --time=96:00:00 \
    --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python39/bin/python -c \"from aces.imaging.mosaic_12m import make_giant_mosaic_cube_sio21; make_giant_mosaic_cube_sio21(channels='all', skip_channel_mosaicing=True, verbose=True,)\"")
