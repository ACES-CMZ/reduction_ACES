"""
parallel_tclean sometimes used `imageconcat` in move mode, but we need everything in `p` (paged, single-file) mode.
"""

import os
import shutil

from casatools import image

from aces import conf

ia = image()

suffixes = "image,weight,image.pbcor,mask,model,pb,psf,residual,sumwt".split(",")

grouppath = 'data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/'

for root, dirs, files in os.walk(os.path.join(conf.basepath, grouppath)):
    for suffix in suffixes:
        if root.endswith(suffix) and any((dirname.endswith(suffix) for dirname in dirs)):
            toconcat = sorted([dirname for dirname in dirs if dirname.endswith(suffix)])
            print(f"Directory {root} appears to be a .move directory, containing dirs {toconcat}")
            ia.imageconcat(outfile=root + ".p",
                           infiles=[os.path.join(root, dirname)
                                    for dirname in toconcat],
                           mode='p')
            shutil.move(root, root + ".move")
            shutil.move(root + ".p", root)
