#!/bin/bash

#argument 1 - The starting compressed folder that includes the db.db for the aligned shapedb database
#argument 2 - Ligand for a shapedb library to be compared against, must be aligned with the ligands/conformers in the database

#get the absolute path to make_db_and_condense_params.sh for use in calling some other scripts in the github repo
# Save the current working directory

original_dir=$(pwd)

# Get the directory where the script is located
script_dir=$(dirname "$(realpath "$0")")

# Change to the script's directory
cd "$script_dir"

# Find the root of the GitHub repository
repo_root=$(git rev-parse --show-toplevel)

# Return to the original directory
cd "$original_dir"


# Now `repo_root` holds the path to the root of the GitHub repository
echo "The root of the GitHub repository is: $repo_root"

#derive the base directory name without the .tar.gz
dirname=$(basename "$1" .tar.gz)
echo $dirname
#untar the ligand confs directory that contains the conformers in sdf and params format, as well as the text file list containing the name of every conformer (lig name + conformer id number) to be used for searching
tar -xzf $1
#move into the directory
cd "$dirname"
#run python script to determine the number of lines in the list file (equal to the number of conformers) and then run the shapedb nnsearch

python $repo_root/shapedb/run_nnsearch.py ${dirname}_lig_name_list.txt $2

#move the resulting text file up and then compress it
#dir_and_sub + "_scored_confs_against_" + sdf_base + ".txt"

# Search for multiple files containing "_scored_confs_against_"
for file in $(find $original_dir/$dirname -type f -name "*_scored_confs_against_*"); do
    echo "Processing file: $file"
    # Do something with $file
	
    mv $file ..

	cd ..

	filename=$(basename "$file")

	#suvo_NN_$(sub_num)_$(directory).tar.gz
	tar -cvf $filename.tar.gz $filename

	rm  $filename

    cd "$dirname"
done

