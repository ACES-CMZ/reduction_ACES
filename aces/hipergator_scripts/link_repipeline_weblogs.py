import os
import glob
from astropy import log

from .. import conf
basepath = conf.basepath

def main():
    for dirname in glob.glob(f"{basepath}/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member*"):
        member = os.path.basename(dirname)

        webpath = f"{basepath}/weblogs-reimaging/{member}/"
        if not os.path.exists(webpath):
            os.mkdir(webpath)

        for pipeline in glob.glob(os.path.join(dirname, "calibrated/working/pipeline*")):
            outpath = os.path.join(webpath, os.path.basename(pipeline))
            if os.path.isdir(pipeline) and not os.path.exists(outpath):
                os.symlink(pipeline, outpath)
            log.info(f'Linked {pipeline} to {outpath}')
