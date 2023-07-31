import os
from tqdm import tqdm
from pathlib import Path
from feather_funcs import *
from cube_utils import *
from reproject_mosaic_funcs import *

"""
- Define the relevant paths for the ACES pipeline
- Choose the line(s) to be processed
- Choose the desired velocity range and resolution
- Choose whether to include the 12M data (if False, only the TP+7M data will be processed)
- Choose the pixel rebinning factor
"""
ACES_ROOTDIR = Path(os.getenv('ACES_ROOTDIR'))
ACES_WORKDIR = Path(os.getenv('ACES_WORKDIR'))
ACES_DATA = Path(os.getenv('ACES_DATA'))

LINE_TABLE = (ACES_ROOTDIR / 'aces/data/tables/linelist.csv')
LINES = ['hnco43']

START_VELOCITY = -100  # km/s
END_VELOCITY = 100  # km/s
VEL_RES = 1  # km/s

INCLUDE_12M = True
REBIN_FACTOR = 2

"""
The following loop will process each line in the LINES list. For each line, the following steps will be performed:
- Read the line table
- Create the feathered cubes
- Crop the image and weight cubes to the desired velocity range and resolution (this step is slow)
- Create a weighted mosaic that is smoothed to the largest beam size in the individual cubes
- Convert the mosaic to K
- Rebin the mosaic (if desired)
"""

for LINE in tqdm(LINES, desc='LINES'):
    LINE_SPWS = read_table(LINE_TABLE)
    create_feathercubes(ACES_WORKDIR, ACES_DATA, ACES_ROOTDIR, LINE_SPWS, LINE, process_12M=INCLUDE_12M)
    crop_cubes(ACES_WORKDIR, START_VELOCITY, END_VELOCITY, VEL_RES, LINE_SPWS, LINE, process_12M=INCLUDE_12M)
    create_weighted_mosaic(ACES_WORKDIR, START_VELOCITY, END_VELOCITY, VEL_RES, LINE, process_12M=INCLUDE_12M)
    cubeconvert_K_kms(ACES_WORKDIR, LINE, START_VELOCITY, END_VELOCITY, VEL_RES, process_12M=INCLUDE_12M)
    rebin(ACES_WORKDIR, LINE, START_VELOCITY, END_VELOCITY, VEL_RES, REBIN_FACTOR, process_12M=INCLUDE_12M)