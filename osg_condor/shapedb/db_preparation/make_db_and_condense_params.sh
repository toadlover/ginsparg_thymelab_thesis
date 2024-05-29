#!/bin/bash

#set -e
#$3 is sub_directory number (0-9)
#$4 is the directory (00000-53084)

#set up pip and import relevant packages (may need to incorporate this into the base container)
#apt-get update
#apt-get install python-pip -y
#pip install numpy
#pip install openbabel


#untar the split_new_named directory that contains the conformers in sdf and params format, as well as the text file list containing the name of every conformer (lig name + conformer id number) to be used for searching
tar -xzvf split_new_named_$3.tar.gz
#move into the directory
cd split_new_named_$3

#delete the single sdf file folder (we don't need it)
rm -drf single_conf_sdfs

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
python3 ../../condense_conformer_params_data.py
#text files of condensed params should be generated. We can delete the .params to save on a lot of space
rm *.params

#back up so we can work on making the database
cd ..

#==========================================================================
#use shapedb to make database (db.db)
#run align.py on the confs_named.sdf file, name as aligned.sdf
python ../align.py confs_named.sdf aligned.sdf

#run shapedb on aligned.sdf to create the database
/pharmit/src/build/shapedb -h -Create -in aligned.sdf -db db.db

#we no longer need confs_named or aligned.sdf
rm confs_named.sdf aligned.sdf

#rename the conformer name text file
mv split_new_named_$3_lig_name_list.txt $4_$3_lig_name_list.txt

#copy the file up a level
#cp $4_$3_lig_name_list.txt ..

#move up a level
cd ..

#rename and compress the folder
mv split_new_named_$3 condensed_params_and_db_$3
tar -czvf condensed_params_and_db_$3.tar.gz condensed_params_and_db_$3

#compress the name list
#tar -czvf $4_$3_lig_name_list.txt.tar.gz $4_$3_lig_name_list.txt

#move both new compressed files into a new folder named after the subjob, which will be sent to LTS
#mkdir $3
#mv $4_$3_lig_name_list.txt.tar.gz condensed_params_and_db_$3.tar.gz $3
