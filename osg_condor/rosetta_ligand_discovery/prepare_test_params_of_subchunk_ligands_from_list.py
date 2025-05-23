#the purpose of this script is to take a list of ligands that correspond to a superchunk and extract the params of all conformers for ligands in a specified subchunk
#all conformers will be put in an uncompressed test_params directory
#note, run this at the same level as the subchunk data (i.e. condensed_params_and_db_#)
#the extract_single_param_from_condensed_file.py script needs to be at this level as well (this will be the case inside of condor jobs)

import os,sys

#input 1: list csv file, must at the very least have the chunk at index 0, subchunk at index 1, and ligand name at index 3
list_file = sys.argv[1]

#input 2: the working chunk (number from 00000-53084, must include leading zeroes up to a 5 digit number)
working_chunk = str(sys.argv[2])

#input 3: the working subchunk (number from 0-9)
working_subchunk = str(sys.argv[3])

#make a test_params folder
os.system("mkdir test_params")

#prepare test_params files
os.system("echo \"## the atom_type_set and mm-atom_type_set to be used for the subsequent paramet\" > test_params/residue_types.txt")
os.system("echo \"TYPE_SET_MODE full_atom\" >> test_params/residue_types.txt")
os.system("echo \"ATOM_TYPE_SET fa_standard\" >> test_params/residue_types.txt")
os.system("echo \"ELEMENT_SET default\" >> test_params/residue_types.txt")
os.system("echo \"MM_ATOM_TYPE_SET fa_standard\" >> test_params/residue_types.txt")
os.system("echo \"ORBITAL_TYPE_SET fa_standard\" >> test_params/residue_types.txt")
os.system("echo \"## Test params files\" >> test_params/residue_types.txt")
#will write params files after

#write additional needed files
os.system("touch test_params/exclude_pdb_component_list.txt")
os.system("touch test_params/patches.txt")

#now, read through the list file and extract all conformer params for all ligands in the given chunk and subchunk
read_file = open(list_file,"r")

for line in read_file.readlines():
	#extract the chunk, subchunk, and ligand name
	chunk = line.split(",")[0].strip()
	subchunk = line.split(",")[1].strip()
	ligand = line.split(",")[2].strip()

	#continue of the chunk or subchunk are a mismatch, which will be most
	if chunk != working_chunk or subchunk != working_subchunk:
		continue

	num_confs = 0

	#read the corresponding condensed params file to identify how many conformers exist for this ligand (up to 15), ideally do not want to overflow
	lig_file = open("condensed_params_and_db_" + working_subchunk + "/single_conf_params/" + ligand + "_shorthand_params.txt", "r")

	for line2 in lig_file.readlines():
		#we can derive the conformer count from the number of commas in the line starting with "_a:", which defines all conformer names
		if line2.startswith("_a:"):
			#this is technically the number + 1, but we need that for iterating over a range anyway
			num_confs = len(line2.split(","))

	#if we are here, then the ligand is for this chunk and subchunk
	#iterate from 1-15 to run the extract_single_param_from_condensed_file.py script to extract all conformers of this ligand
	for i in range(1,num_confs):
		os.system("python extract_single_param_from_condensed_file.py condensed_params_and_db_" + working_subchunk + "/single_conf_params/" + ligand + "_shorthand_params.txt " + str(i) + " " + ligand + "_" + str(i))

		#write the newly made params file to the residue_types file and then move it to test_params
		os.system("ls *.params >> test_params/residue_types.txt")
		os.system("mv *.params test_params")
		os.system("rm *.params")

	#we are now done with the ligand and can move to the next