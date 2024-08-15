#!/bin/bash

#set -e
#$1 is the file

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
dirname=$(basename "$tarfile" .tar.gz)

#untar the ligand confs directory that contains the conformers in sdf and params format, as well as the text file list containing the name of every conformer (lig name + conformer id number) to be used for searching
tar -xzf $1
#move into the directory
cd $dirname

#delete the single sdf file folder (we don't need it)
rm -drf single_conf_sdfs

#prepare anaconda
whoami
. /opt/miniconda3/bin/activate
env "PATH=$PATH" conda update conda
/opt/miniconda/bin/conda init bash
. ~/.bashrc
/opt/miniconda/bin/conda activate
/opt/miniconda/bin/conda list

#==========================================================================
#params handling
#move into the params directory and then run the params condense script
cd single_conf_params
python3 $repo_root/params_file_compression/condense_conformer_params_data.py
#text files of condensed params should be generated. We can delete the .params to save on a lot of space
rm *.params

#back up so we can work on making the database
cd ..

#==========================================================================
#use shapedb to make database (db.db)
#run align.py on the confs_named.sdf file, name as aligned.sdf
python $repo_root/shapedb/align.py confs_named.sdf aligned.sdf

#run shapedb on aligned.sdf to create the database
/pharmit/src/build/shapedb -h -Create -in aligned.sdf -db db.db

#we no longer need confs_named or aligned.sdf
rm confs_named.sdf aligned.sdf

#move up a level
cd ..

#recompress the folder

tar -czf $1 $dirname

#compress the name list
#tar -czvf $4_$3_lig_name_list.txt.tar.gz $4_$3_lig_name_list.txt

#move both new compressed files into a new folder named after the subjob, which will be sent to LTS
#mkdir $3
#mv $4_$3_lig_name_list.txt.tar.gz condensed_params_and_db_$3.tar.gz $3
