import os
import re
import sys
import glob
import json
import pandas as pd

"""
Crude script to pull from default_tclean_commands.json file to grab all region
names and corresponding 12m & 7m MOUS IDs. This information is written out to a
CSV file containing the following columns:
- Region name (with unnecessary text trimmed, e.g. "_updated"/"_03_TM1")
- Original, untrimmed 12m region name (for querying .json file in cleaning script)
- 12m MOUS ID (e.g. uid___A001_X1590_X30aa)
- Original, untrimmed 7m region name
- 7m MOUS ID
"""

rootdir = os.getenv('ACES_ROOTDIR')
workdir = os.getenv('ACES_WORKDIR')

with open(f"{rootdir}/aces/pipeline_scripts/default_tclean_commands.json", "r") as fh:
    df = pd.read_json(fh)

# Grab all region names and MOUS IDs from .json file
id_12m = []
reg12m = []
id_7m  = []  # noqa: E221
reg7m  = []  # noqa: E221
for col in df.columns:
    if 'TM1' in col:
        reg12m.append(col)
        id_12m.append(df[col]['tclean_cube_pars']['spw25']['imagename'].removesuffix('.s38_0.Sgr_A_star_sci.spw25.cube.I.iter1'))
    if '7M' in col:
        reg7m.append(col)
        id_7m.append(df[col]['tclean_cube_pars']['spw16']['imagename'].removesuffix('.s38_0.Sgr_A_star_sci.spw16.cube.I.iter1'))

# Create pandas dataframes from above lists
sb_names_12m = pd.DataFrame(
    {'Region': reg12m,
     'Reg_original_12m': reg12m,
     'twelve_m_ID': id_12m,
    })  # noqa: E124

sb_names_7m = pd.DataFrame(
    {'Region': reg7m,
     'Reg_original_7m': reg7m,
     'seven_m_ID': id_7m
    })  # noqa: E124

# Discard repeated regions (i.e. those with wrong freq shift)
updated_12m = [s for s in sb_names_12m['Region'] if "updated" in s]
for i in range(len(updated_12m)):
    updated_12m[i] = re.sub('_updated', '', updated_12m[i])
sb_12m = sb_names_12m[~sb_names_12m['Region'].isin(updated_12m)]

updated_7m = [s for s in sb_names_7m['Region'] if "updated" in s]
for i in range(len(updated_7m)):
    updated_7m[i] = re.sub('_updated', '', updated_7m[i])
sb_7m = sb_names_7m[~sb_names_7m['Region'].isin(updated_7m)]

# Tidy up region names, merge DFs on region name, and export as CSV
sb_12m['Region'] = sb_12m['Region'].str.replace('_updated', '')
sb_12m['Region'] = sb_12m['Region'].str.replace('_03_TM1', '')
sb_7m['Region']  = sb_7m['Region'].str.replace('_updated', '')  # noqa: E221
sb_7m['Region']  = sb_7m['Region'].str.replace('_03_7M', '')  # noqa: E221
sb_names         = pd.merge(sb_12m, sb_7m, on='Region', how='outer').dropna().sort_values('Region')  # noqa: E221
sb_names.to_csv(workdir+'/'+'aces_SB_names.csv', index=False)  # noqa: E226
