#10-9-24
#The purpose of this script is to take in a text document containing a list of pdb files (one file per line, paths can be included), extract the ligand and append all ligands in all files into a mol2 file.
#A csv file containing the ligand names and smiles strings will also be created
#the input text file, desired name for output sdf file, and desired name for output csv file are all mandatory arguments for this script
#this script requires an environment with rdkit (for SMILES) and openbabel (for pdb->sdf conversion)

#this script is expected to be used on placed ligand files from the Rosetta LigandDiscoverySearch
#name of input pdb file in the text list should look something like: 4s0v_receptor_only_Z3360656891_97_0.pdb where the ligand name is between the 3rd and 4th to last underscores.
#this is an issue if your ligand name has an underscore in it, but Enamine and Zinc ligand named do not have underscores

#example run command: python prepare_ligands_from_pdb_for_enamine_order.py antagonist_drugs_file_list.txt antagonists_to_order.sdf antagonists_to_order.csv

import os,sys
#use Chem in rdkit to derive SMILES string of placed ligands
from rdkit import Chem

#take input arguments
input_text_file_name = sys.argv[1]
output_sdf_file_name = sys.argv[2]
output_csv_file_name = sys.argv[3]

#create file streams for the input file and output csv file
input_file = open(input_text_file_name,"r")
output_csv = open(output_csv_file_name,"w")

#run through the input file
for line in input_file.readlines():
	#strip any newline off the end of the line
	stripped_newline = line.strip("\n")

	#derive the ligand name
	lig_name = stripped_newline.split("_")[len(stripped_newline.split("_")) - 3]

	print("On ligand: " + lig_name)

	#extract the ligand as HETATM lines to a temporary pdb file to be converted to sdf
	os.system("grep HETATM " + stripped_newline + " > temp.pdb")

	#use openbabel to convert temp.pdb to temp.sdf
	os.system("obabel temp.pdb -O temp.sdf")

	#change the name in temp.sdf to be the ligand name using sed (it will previously be temp.pdb)
	os.system("sed -i 's/temp.pdb/" + lig_name + "/' temp.sdf")

	#append temp.sdf to the output sdf file
	os.system("cat temp.sdf >> " + output_sdf_file_name)

	#now derive the smiles string of the ligand

	smiles = ""

	#derive the smiles string
	supplier = Chem.SDMolSupplier("temp.sdf")
	for mol in supplier:
		if mol is not None:  # Check for valid molecule
			smiles = Chem.MolToSmiles(mol)

	#write the smiles to the csv file with the ligand name
	output_csv.write(lig_name + "," + smiles + "\n")

#cleanup

#remove temp.pdb and temp.mol2
os.system("rm temp.pdb temp.sdf")

output_csv.close()
input_file.close()