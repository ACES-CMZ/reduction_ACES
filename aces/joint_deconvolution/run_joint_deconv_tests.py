import os
from tqdm import tqdm
from pathlib import Path
from joint_deconv_funcs import do_joint_deconvolution, read_table, process_string

ACES_ROOTDIR = Path(os.getenv('ACES_ROOTDIR'))
ACES_WORKDIR = Path(os.getenv('ACES_WORKDIR'))
ACES_DATA = Path(os.getenv('ACES_DATA'))
LINE_TABLE = (ACES_ROOTDIR / 'aces/data/tables/linelist.csv')

# Set up the parameters for joint deconvolution
LINES = ['hnco43']
V_START = '-20 km/s'
V_WIDTH = '1 km/s'
NCHAN = 120
REGION = ['ao']

DEEP_CLEAN = True
RELAXED_MASKING = True
TP_STARTMODEL = False
#####################################################


for LINE in tqdm(LINES, desc='LINES'):
    LINE_SPWS = read_table(LINE_TABLE)
    RESTFREQ = LINE_SPWS[process_string(LINE)]['restfreq']
    do_joint_deconvolution(ACES_WORKDIR, ACES_DATA, ACES_ROOTDIR, REGION, LINE_SPWS, LINE, RESTFREQ, V_START, V_WIDTH, NCHAN, DEEP_CLEAN, RELAXED_MASKING, TP_STARTMODEL)