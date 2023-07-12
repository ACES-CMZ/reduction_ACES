#!/bin/bash
# This is an example template for using the Globus CLI to transfer data from HPG to a local machine
# This example is for transferring the TP fits files, and uses the aces_SB_uids.csv file to get the list of SB uids
# This is designed to work with a directory structure that directly mirrors the HPG directory structure
# Run by chmod +x globus_cli.sh
# Run by ./globus_cli.sh

# Define endpoint IDs and target paths - HPG
# hpg=insert_hpg_endpoint_id_here
hpg=8851ceea-6ae4-4c1d-8407-ac7a70ee0274 
hpg_pth=/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/

# Define endpoint IDs and target paths - LOCAL
# local=insert_local_endpoint_id_here
local=56bc2744-f0ce-11ed-9a77-83ef71fbf0ae
local_pth=/export/data1/abarnes/alma-aces/feather/

# UID csv file
csv_file='./../../tables/aces_SB_uids.csv'
column_number=4

# Remove any existing file named aces_fits.txt
rm aces_fits.txt

echo "[INFO] Adding transfer files to list"
# Read the CSV file and extract the specified column values
# Skip the header row and print the value of the specified column for each subsequent row
awk -F ',' -v col=$column_number '{
    if (NR > 1) {  # Skip the header row
        print $col  # Print the value of the specified column
    }
}' "$csv_file" | while read -r uid; do

    # Define the source path by concatenating the HPG path and the UID
    src_path=$hpg_pth"member.uid___A001_"$uid"/product/"

    # Print the command line with the filename
    command="$hpg:~$src_path"
    echo "[INFO] Adding file: $command"

    # Append the paths of the files to be transferred to the aces_fits.txt file
    # Filter the files using the specified pattern (e.g. ~*spw31.cube.I.pbcor.fits*)
    # Use xargs to pass each file to the echo command and append the path to the text file
    globus ls "$command" --filter '~*cube.I.sd.fits' | xargs -i echo $src_path"{}" >> aces_fits.txt
done

echo "[INFO] Executing transfer commands"
# Transfer the files
while IFS= read -r fn;
do
  if [ ! -f "~/$fn" ]; then
    transfer_command="globus transfer \"$hpg:~$fn\" \"$local:$local_pth$fn\""
    echo "[INFO] Executing transfer command: $transfer_command"
    eval "$transfer_command"
  fi
done < aces_fits.txt