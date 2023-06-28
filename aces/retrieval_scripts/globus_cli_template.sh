#!/bin/bash
# This is an example template for using the Globus CLI to transfer data from HPG to a local machine
# This example is for transferring the 12m HNCO fits files, and uses the aces_SB_uids.csv file to get the list of SB uids
# This is designed to work with a directory structure that directly mirrors the HPG directory structure

# Define endpoint IDs and target paths
hpg=insert_hpg_endpoint_id_here
hpg_pth=/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/
local=insert_local_endpoint_id_here
local_pth=~/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/

csv_file='./../data/tables/aces_SB_uids.csv'
column_number=2

rm aces_12m_HNCO_fits.txt
awk -F ',' -v col=$column_number '{
    if (NR > 1) {  # Skip the header row
        print $col  # Print the value of the specified column
    }
}' "$csv_file" | while read -r uid; do
    # Define the source path
    src_path=$hpg_pth"member.uid___A001_"$uid"/"
    # Define the destination path
    dst_path=$local_pth"member.uid___A001_"$uid"/"
    mkdir -p "$dst_path"
    mkdir -p $dst_path"product/"
    # Write the paths of the files to be transferred to a text file
    globus ls $hpg:\~$src_path"product/" --filter '~*spw31.cube.I.pbcor.fits*' | xargs -i echo $src_path"product/{}" >> aces_12m_HNCO_fits.txt
done

# Transfer the files
while IFS= read -r fn;
do
  if [ ! -f "~/$fn" ]; then
    globus transfer "$hpg:~$fn" "$local:~/$fn"
  fi
done < aces_12m_HNCO_fits.txt