#!/bin/bash

#$3 is the superchunk file name
#$4 is the chunk number
#$5 is the subchunk number
#$6 is the arg file name



begin_time="$(date +%s.%N)"
echo $begin_time

echo $3 $4 $5 $6

#unzip the condensed_params_and_db_#.tar.gz
tar -xzf condensed_params_and_db_$5.tar.gz

echo "post-decompress ls"
ls

chmod -R 777 *

echo "extracting all conformer params for specified ligands"

python prepare_test_params_of_subchunk_ligands_from_list.py $3 $4 $5

ls

echo "File count in unzipped test params"
ls test_params | wc -l

cat test_params/residue_types.txt

setup_operations_end="$(date +%s.%N)"
echo $setup_operations_end

echo running rosetta

#/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @$3 > rosetta_output.txt
/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @$6

rosetta_end="$(date +%s.%N)"
echo $rosetta_end

ls

mkdir placements
#mv *delta*pdb placements
rm minipose.pdb
mv *pdb placements

#compress pdb folder
tar -czf placements.tar.gz placements

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
