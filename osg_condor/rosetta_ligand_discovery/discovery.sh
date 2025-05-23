#!/bin/bash

#$3 is the ligand name
#$4 is the conformer number (0-14)
#$5 is the source chunk
#$6 is the source subchunk
#$7 is the ligand used by shapedb
#$8 is the shapedb score
#$9 is the custom arg file name inputted in the .sub file

begin_time="$(date +%s.%N)"
echo $begin_time

chmod -R 777 *

echo $3 $4 $5 $6 $7 $8

echo Initial ls
ls

echo Making test_params

#make test_params directory
mkdir test_params
#make residue_types file
echo "## the atom_type_set and mm-atom_type_set to be used for the subsequent paramet" > residue_types.txt
echo "TYPE_SET_MODE full_atom" >> residue_types.txt
echo "ATOM_TYPE_SET fa_standard" >> residue_types.txt
echo "ELEMENT_SET default" >> residue_types.txt
echo "MM_ATOM_TYPE_SET fa_standard" >> residue_types.txt
echo "ORBITAL_TYPE_SET fa_standard" >> residue_types.txt
echo "## Test params files" >> residue_types.txt
echo "$3.params" >> residue_types.txt
mv residue_types.txt test_params

touch test_params/exclude_pdb_component_list.txt
touch test_params/patches.txt

echo Unzipping tar
#untar the directory that contains the db.db so that the nnsearch can be run
tar -xzvf condensed_params_and_db_$6.tar.gz
#move into the directory
pwd
ls
cd condensed_params_and_db_$6/single_conf_params
pwd
ls

echo running python script
which python
#run python script to decompress the corresponding param file for the right ligand and conformer
python ../../extract_single_param_from_condensed_file.py $3_shorthand_params.txt $4 $3

echo test print for confirming params generation
cat *params
ls *params

##TEMPORARY FIX BECAUSE CONDENSED PARAMS FILES CAN'T BE READ BY ROSETTA
#sed "s/ /  /g" $3.params > temp.params
#mv temp.params $3.params

python ../../fix_condensed_param_file_spacing.py $3.params
mv fixed_$3.params  $3.params

cat $3.params

#move params up 2 levels and move yourself up
mv $3.params ../../test_params



cd ../..

setup_operations_end="$(date +%s.%N)"
echo $setup_operations_end

echo running rosetta
#run rosetta ligand discovery on all 5 anchor residues
#/main/source/bin/ligand_discovery_search_protocol.linuxgccrelease @argsASN227
#/main/source/bin/ligand_discovery_search_protocol.linuxgccrelease @argsGLN86
#/main/source/bin/ligand_discovery_search_protocol.linuxgccrelease @argsHIS253
#/main/source/bin/ligand_discovery_search_protocol.linuxgccrelease @argsTHR63
#/main/source/bin/ligand_discovery_search_protocol.linuxgccrelease @argsTYR257

/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @$9 > rosetta_output.txt

rosetta_end="$(date +%s.%N)"
echo $rosetta_end

#move output pdbs to location for zipping and storing
mkdir $3_$4_$7
mv *.pdb $3_$4_$7

#python initial_score_breakdown.py $8
python initial_score_breakdown_from_rosetta_output.py $8 $3 $4 $5 $6

#compress pdb folder
tar -czvf $3_$4_$7.tar.gz $3_$4_$7

post_operations_end="$(date +%s.%N)"
echo $post_operations_end

#make file with the times

#can't use bc because the rosetta container lacks bc. Going to try using awk operations
#startup operations
#startup_operations=$(echo "$setup_operations_end-$begin_time" |bc)
startup_operations=$(echo $setup_operations_end $begin_time | awk '{print $1 - $2}')
#rosetta call
#rosetta_operations=$(echo "$rosetta_end-$setup_operations_end" |bc)
rosetta_operations=$(echo $rosetta_end $setup_operations_end | awk '{print $1 - $2}')
#post operations
#end_operations=$(echo "$post_operations_end-$rosetta_end" |bc)
end_operations=$(echo $post_operations_end $rosetta_end | awk '{print $1 - $2}')

echo $startup_operations
echo $rosetta_operations
echo $end_operations

echo startup:$startup_operations,rosetta:$rosetta_operations,end:$end_operations > job_times.txt
