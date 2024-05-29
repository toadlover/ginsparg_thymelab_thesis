#!/bin/bash

#set -e
#$3 is sub_directory number (0-9)
#$4 is the directory (00000-53084)

#untar the directory that contains the db.db so that the nnsearch can be run
tar -xzvf $3.tar.gz
#move into the directory
pwd
ls
cd condensed_params_and_db_0
pwd
ls
#run python script to determine the number of lines in the list file (equal to the number of conformers) and then run the shapedb nnsearch

python ../run_nnsearch.py 00000_0_lig_name_list.txt ../suvo_shifted.sdf

#move the resulting text file up and then compress it
#dir_and_sub + "_scored_confs_against_" + sdf_base + ".txt"
mv 00000_0_scored_confs_against_suvo_shifted.txt ..
cd ..

#suvo_NN_$(sub_num)_$(directory).tar.gz
tar -czvf suvo_NN_$3.tar.gz 00000_0_scored_confs_against_suvo_shifted.txt

