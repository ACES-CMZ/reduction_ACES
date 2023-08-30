

incr=50

for ch in `seq 0 ${incr} 400`; do
    ch1=$ch
    ch2=$((${ch1} + ${incr}))
    sbatch --job-name=aces_sio21_mosaic_ch${ch1}to${ch2} \
        --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_sio21_mosaic_ch${ch1}to${ch2}_%j.log  \
        --account=astronomy-dept --qos=astronomy-dept-b \
        --ntasks=8 --nodes=1 --mem=32gb --time=96:00:00 \
        --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python39/bin/python -c \"from aces.imaging.mosaic_12m import sio21_cube_mosaicing; sio21_cube_mosaicing(channels=range(${ch1},${ch2}), verbose=True)\""
done
