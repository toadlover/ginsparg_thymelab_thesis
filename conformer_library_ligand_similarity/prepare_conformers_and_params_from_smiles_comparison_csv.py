#the purpose of this script is to take in a csv file that has the smiles strings of ligands of interest
#we take the smiles, derive sdf files of the smiles, and use conformator to build a diverse conformer set (use default cutoff of 150), and make rosetta params of the conformers to investigate with Rosetta
#each ligand will get its own folder, and will be uploaded to a bucket location of choice

#csv line format:
#chunk,subchunk,ligname,similarityscore,smiles

#needed inputs
#results csv file
#path to conformator executable
#path to working copy of Rosetta molfile_to_params.py
#bucket location to push compressed test_params directories

#package imports
#this utilizes rdkit, as well as having conformator, Rosetta, and s3cmd
#actually, I don't think we need rdkit
import os,sys
#from rdkit import Chem
#from rdkit.Chem import AllChem

#read in mandatory inputs in order

#working location
working_location = sys.argv[1]

#results csv file, can include path
csv_file = sys.argv[2]

#conformator executable
#path to and conformator executable, i.e. /pi/summer.thyme-umw/2024_intern_lab_space/conformator_1.2.1/conformator
conformator_executable = sys.argv[3]

#molfile_to_params executable/python script path and name
m_to_p_executable = sys.argv[4]

#move to the working location
os.chdir(working_location)

#open the csv file and begin to process the data
read_file = open(csv_file,"r")

#iterate through each line in the file and process each respective ligand
for line in read_file.readlines():
	#extract the ligand name and smiles string, and write a smiles file
	lig_name = line.split(",")[2]
	smiles = line.split(",")[4].strip()

	smiles_file = open("temp.smi", "w")
	smiles_file.write(smiles)
	smiles_file.close()

	#run conformator on the temp smiles file
	os.system(conformator_executable + " -i temp.smi" + " -o " + lig_name + "_confs.sdf --keep3d --hydrogens -v 0")

	#now make a directory for the conformers
	os.system("mkdir " + lig_name)

	#move the new conformers file into the directory
	os.system("mv " + lig_name + "_confs.sdf " + lig_name)

	#enter the directory and process this file
	os.chdir(lig_name)

	#use obabel from command line to split the file
	os.system("obabel " + lig_name + "_confs.sdf -O " + lig_name + "_.sdf -m")

	#delete the original
	os.system("rm " + lig_name + "_confs.sdf ")

	#run through each file and fix the name line, since it would be blank by default
	#run through each newly made file and adjust the ligand name in the file
	for r2,d2,f2 in os.walk(os.getcwd()):
		for single_file in f2:
			#read the file and write to a temporary copy
			read_file = open(single_file,"r")
			write_file = open("temp.sdf", "w")

			#line counter, we are only interested in line 1
			line_counter = 0

			for line2 in read_file.readlines():
				if line_counter == 0:
					write_file.write(single_file.split(".")[0] + "\n")
				else:
					write_file.write(line2)

				line_counter = line_counter + 1

			read_file.close()
			write_file.close()

			#write the temp over the original
			os.system("mv temp.sdf " + single_file)

			#make a params file
			os.system("python " + m_to_p_executable + " " + single_file + " -n " + single_file.split(".")[0] + " --keep-names --long-names --clobber --no-pdb")

	#now, remove the sdf files and make a compressed test_params folder
	os.system("rm *sdf")
	os.system("mkdir test_params")
	os.system("mv *params test_params")

	#move into the test_params folder and set it up
	os.chdir("test_params")

	#touch necessary files
	os.system("touch exclude_pdb_component_list.txt patches.txt")

	#make the residue_types file
	res_types_file = open("residue_types.txt", "w")
	res_types_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
	res_types_file.write("TYPE_SET_MODE full_atom\n")
	res_types_file.write("ATOM_TYPE_SET fa_standard\n")
	res_types_file.write("ELEMENT_SET default\n")
	res_types_file.write("MM_ATOM_TYPE_SET fa_standard\n")
	res_types_file.write("ORBITAL_TYPE_SET fa_standard\n")
	res_types_file.write("## Test params files\n")
	#write all params file names to this file
	for r,d,f in os.walk(os.getcwd()):
		for file in f:
			res_types_file.write(file + "\n")

	#go back up
	os.chdir("..")

	#compress the directory and toss the original
	os.system("tar -czf  test_params.tar.gz test_params")

	#end behavior, go back up a directory
	os.chdir("..")


	#53079,9,PV-006710095263,1.0,CC1CN(C(=O)CCc2ccc3ccccc3c2O)CCN1C(=O)c1cc[nH]n1
	#
#python /pi/summer.thyme-umw/2024_intern_lab_space/ari_work/ginsparg_thymelab_thesis/conformer_library_ligand_similarity/prepare_conformers_and_params_from_smiles_comparison_csv.py /pi/summer.thyme-umw/2024_intern_lab_space/ari_work/10k_drug_27_ligands/test 0_drug_27_best_sorted_test_chiral.csv /pi/summer.thyme-umw/2024_intern_lab_space/conformator_1.2.1/conformator /pi/summer.thyme-umw/2024_intern_lab_space/rosetta/source/scripts/python/public/molfile_to_params.py