#this is a modified version of the make conformer and params files from a single csv file that contains smiles data
#the file must be formatted with the smiles string in the first comma index, and the ligand name in the second. extra data can be stored beyond, but will not be interacted with
#this will also make a functional test_params directory containing all params files made at the end

import os,sys

#csv file with smiles strings
starting_file = sys.argv[1]

#license key used to activate conformator
license_key = sys.argv[2]

#split name to remove extension
#if there is a path, remove that too to work in the current location
starting_file_no_ext = starting_file.split("/")[len(starting_file.split("/")) - 1].split(".")[0]

#make a test_params directory to perform all operations in and then move into it
os.system("mkdir test_params" )
os.system("cp " + starting_file + " test_params")
os.chdir(test_params)

#activate conformator
os.system("/conformator_for_container/conformator_1.2.1/conformator --license " + license_key)

#iterate over each line of the csv file to make conformers and params
smiles_csv = open(starting_file_no_ext,"r")
for line in smiles_csv.readlines():
	#extract the smiles string and ligand name
	smiles_string = line.split(",")[0].strip()
	lig_name = line.split(",")[1].strip()

	#write the smiles string to a smiles file
	os.system("echo \"" + smiles_string + "\" > smiles.smi")

	#run conformator on the starting file
	os.system("/conformator_for_container/conformator_1.2.1/conformator -i smiles.smi -o " + lig_name + ".sdf --keep3d --hydrogens -n 15 -v 0")

	#use babel to break the confs.sdf file into individual sdf files; also apply unique names to the ligand name (per conformer)
	os.system("babel " + lig_name + ".sdf " + lig_name + ".sdf -m")

	#iterate over each conformer sdf and make params for it
	for r,d,f in os.walk(os.getcwd()):
		for file in if:
			#if it is a ligand conformer sdf
			if lig_name in file and file.endswith(".sdf"):
				#run molfile to params
				#make a params file of the unique file
				os.system("python /conformator_for_container/molfile_to_params.py " + file + " -n " + file.split(".")[0] + " --keep-names --long-names --clobber --no-pdb")

	#with all params made, clean up by removing all sdf and smi files
	os.system("rm *smi *sdf")

#make the other necessary test_params files
os.system("touch exclude_pdb_component_list.txt")
os.system("touch patches.txt")

#residue types file
os.system("echo \"## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\" >> residue_types.txt")
os.system("echo \"TYPE_SET_MODE full_atom\" >> residue_types.txt")
os.system("echo \"ATOM_TYPE_SET fa_standard\" >> residue_types.txt")
os.system("echo \"ELEMENT_SET default\" >> residue_types.txt")
os.system("echo \"MM_ATOM_TYPE_SET fa_standard\" >> residue_types.txt")
os.system("echo \"ORBITAL_TYPE_SET fa_standard\" >> residue_types.txt")
os.system("echo \"## Test params files\" >> residue_types.txt")
#write all the params files to the file
os.system("ls *.params >> residue_types.txt")