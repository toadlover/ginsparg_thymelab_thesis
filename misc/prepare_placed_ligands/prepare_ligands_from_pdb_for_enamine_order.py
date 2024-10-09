#10-9-24
#The purpose of this script is to take in a text document containing a list of pdb files (one file per line, paths can be included), extract the ligand and append all ligands in all files into a mol2 file.
#A csv file containing the ligand names and smiles strings will also be created
#the input text file, desired name for output mol2 file, and desired name for output csv file are all mandatory arguments for this script
#this script requires an environment with rdkit (for SMILES) and openbabel (for pdb->mol2 conversion)

#this script is expected to be used on placed ligand files from the Rosetta LigandDiscoverySearch
#name of input pdb file in the text list should look something like: 4s0v_receptor_only_Z3360656891_97_0.pdb where the ligand name is between the 3rd and 4th to last underscores.
#this is an issue if your ligand name has an underscore in it, but Enamine and Zinc ligand named do not have underscores

#example run command: python prepare_ligands_from_pdb_for_enamine_order.py antagonist_drugs_file_list.txt antagonists_to_order.mol2 antagonists_to_order.csv

import os,sys
#use Chem in rdkit to derive SMILES string of placed ligands
from rdkit import Chem

#take input arguments
input_text_file_name = sys.argv[1]
output_mol2_file_name = sys.argv[2]
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

	print(lig_name)