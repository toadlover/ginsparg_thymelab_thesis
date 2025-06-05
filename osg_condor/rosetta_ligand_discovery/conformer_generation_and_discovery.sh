#!/bin/bash

#$3 is the conformator license string
#$4 is the arg file name
#$5 is the name of the starting pdb file, we want to delete it so it does not exist as memory overhead


begin_time="$(date +%s.%N)"
echo $begin_time

echo $3 $4

#for the ease of running, rename the ligand csv file to ligs.csv; there should only by the single input csv file
mv *.csv ligs.csv

#run /ginsparg_thymelab_thesis/conformator/make_conformers_and_params_from_smiles_file.py to make conformer params from the inputted csv file
python make_conformers_and_params_from_smiles_file.py ligs.csv $3

#cd /

echo "File count in unzipped test params"
ls test_params | wc -l

cat test_params/residue_types.txt

setup_operations_end="$(date +%s.%N)"
echo $setup_operations_end

echo running rosetta

#/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @$3 > rosetta_output.txt
/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @$4

rosetta_end="$(date +%s.%N)"
echo $rosetta_end

ls

mkdir placements
#mv *delta*pdb placements
rm minipose.pdb $5
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
