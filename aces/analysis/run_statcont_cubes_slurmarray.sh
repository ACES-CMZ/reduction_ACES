
basepath=/orange/adamginsburg/ACES/
filenames=${basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_*/calibrated/working/*cube*.image.pbcor.fits
nfiles=$(ls -1 ${basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_*/calibrated/working/*cube*.image.pbcor.fits | wc -l)
echo "Working on ${nfiles} files"

flist=$(/blue/adamginsburg/adamginsburg/miniconda3/envs/python312/bin/python -c "from aces.analysis.statcont_cubes import get_file_numbers; rslt=get_file_numbers(); print(','.join(map(str,rslt)))" | tail -1)
echo "Sub-list is $flist"

if [ $flist ]; then

    jobid=$(sbatch --job-name=aces_statcont_arr \
        --output=/blue/adamginsburg/adamginsburg/ACES/logs/aces_statcont_arr_%j_%A_%a.log  \
        --array=$flist \
        --account=astronomy-dept --qos=astronomy-dept-b \
        --ntasks=16 --nodes=1 --mem=64gb --time=96:00:00 --parsable \
        --wrap "/blue/adamginsburg/adamginsburg/miniconda3/envs/python312/bin/python -c \"from aces.analysis.statcont_cubes import main; main()\"")

    echo "Job started is $jobid"

fi

#for fn in filenames; do
#    if [ ! -f ${fn/image.pbcor.fits/image.pbcor.statcont.cont.fits}]; then
#        export CUBEFILENAME=$fn
#        basefn=${fn%.*}
#
#    fi
#done
