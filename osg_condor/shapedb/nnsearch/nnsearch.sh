#!/bin/bash

#set -e
#$3 is sub_directory number (0-9)
#$4 is the directory (00000-53084)
#$5 is the file prefix
#$6 is the file type (.sdf or .mol2 ideally)

#untar the directory that contains the db.db so that the nnsearch can be run
tar -xzvf condensed_params_and_db_$3.tar.gz
#move into the directory
pwd
ls
cd condensed_params_and_db_$3 
pwd
ls
#run python script to determine the number of lines in the list file (equal to the number of conformers) and then run the shapedb nnsearch

python ../run_nnsearch.py $4_$3_lig_name_list.txt ../$5.$6

#move the resulting text file up and then compress it
#dir_and_sub + "_scored_confs_against_" + sdf_base + ".txt"
mv $4_$3_scored_confs_against_$5.txt ..
cd ..

#suvo_NN_$(sub_num)_$(directory).tar.gz
tar -czvf $5_NN_$3_$4.tar.gz $4_$3_scored_confs_against_$5.txt

rm  $4_$3_scored_confs_against_$5.txt
