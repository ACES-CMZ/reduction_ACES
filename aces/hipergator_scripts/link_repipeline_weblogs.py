import os
import glob

for dirname in glob.glob("/orange/adamginsburg/ACES/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member*"):
    member = os.path.basename(dirname)

    webpath = f"/orange/adamginsburg/web/secure/ACES/weblogs-reimaging/{member}/"
    if not os.path.exists(webpath):
        os.mkdir(webpath)

    for pipeline in glob.glob(os.path.join(dirname, "calibrated/working/pipeline*")):
        outpath = os.path.join(webpath, os.path.basename(pipeline))
        if os.path.isdir(pipeline) and not os.path.exists(outpath):
            os.symlink(pipeline, outpath)
