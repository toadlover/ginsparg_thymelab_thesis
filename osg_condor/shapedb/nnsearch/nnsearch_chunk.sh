#!/bin/bash

#this runs on (up to) all 10 sub chunks in a chunk, instead of just an individual sub chunk

#set -e
#$3 is the number if sub-chunks - 1 (0-9) # no longer subchunnk, but rather a value to be used to provide a range for a for loop; honestly probably don't even need to use this, just don't use the chunk search on the 53084 chunk and everything should work otherwise
#$4 is the directory (00000-53084)
#$5 is the file prefix
#$6 is the file type (.sdf or .mol2 ideally)

#for i in {1..5}; do echo "Welcome ${i} times"; done

#################################################

#for i in {1..5}
#do
#echo "Welcome ${i} times"
#Done

for i in {0..9}
do
	#untar the directory that contains the db.db so that the nnsearch can be run
	tar -xzvf condensed_params_and_db_${i}.tar.gz
	#move into the directory
	pwd
	ls
	cd condensed_params_and_db_${i}
	pwd
	ls
	#run python script to determine the number of lines in the list file (equal to the number of conformers) and then run the shapedb nnsearch

	python ../run_nnsearch.py $4_${i}_lig_name_list.txt ../$5.$6

	#move the resulting text file up and then compress it
	#dir_and_sub + "_scored_confs_against_" + sdf_base + ".txt"
	mv $4_${i}_scored_confs_against_$5.txt ..
	cd ..

	#suvo_NN_$(sub_num)_$(directory).tar.gz
	tar -czvf $5_NN_${i}_$4.tar.gz $4_${i}_scored_confs_against_$5.txt

	rm  $4_${i}_scored_confs_against_$5.txt
done
