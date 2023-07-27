# Import the necessary libraries
import os
import glob
from tqdm import tqdm
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
from casatasks import imhead, exportfits, imtrans, feather, imreframe
from feather_12m_7m_TP_cube_mosaic_modsfeather import *
from feather_12m_7m_TP_cube_mosaic_modsreprojectcasa import *
from feather_12m_7m_TP_cube_mosaic_modssavefits import *
from astropy.table import Table 


######## ---------------------------------- ##### 
# Set the relevant paths
ACES_DATA = Path('/export/data1/abarnes/alma-aces/feather/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/')
ACES_WORKDIR = Path('/export/data1/abarnes/alma-aces/feather/working/')

# Set path to spectral line table
filename = '../../tables/ACESSpectralLineCoverage.csv'

#Select the molecule/SPW based on the dictionary below, and the code will do the rest.
lines = ['cs21','hnco43','h13cn10', 'ch3ch2cn', 'h13co+10', 'hn13c10', 'sio21', 'hco+10', 'hc3n1110', 'so4544', 'so3221']
######## ---------------------------------- ##### 

for line in tqdm(lines, desc = 'Lines'):

	""" 
		Reads the provided CSV file into a dictionary, where the keys are molecular line names and 
		the values are dictionaries containing the corresponding spectral window (spw) values for 12m, 7m, and TP.
	"""
	# line_spws = read_table(filename)

	""" 
		For each observation listed in the csv file, generates combined cubes (data structures representing 3D images)
		for the total power, 7m, and 12m antennas. The cubes are saved in a directory for each observation.
		The directory structure is created within the working directory (ACES_WORKDIR).
		The raw data for the observations is read from the directory specified by ACES_DATA.
	"""
	# create_feathercubes(ACES_WORKDIR, ACES_DATA, line_spws, line)

	""" 
		This function is responsible for creating a mosaic (a large image constructed from smaller images)
		of the total power, 7m, and 12m cubes for each observation.
		The function will calculate the weighted mean of the overlapping regions between the cubes
		and save the resulting mosaic in the working directory (ACES_WORKDIR).
	"""	
	# create_weighted_mosaic(ACES_WORKDIR, line)
	# cubeconvert_K_kms(ACES_WORKDIR, line)
	rebin(ACES_WORKDIR, line)
