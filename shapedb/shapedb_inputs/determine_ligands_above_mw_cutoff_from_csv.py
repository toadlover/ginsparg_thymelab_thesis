#the purpose of this script is to look at a csv file that results from a smiles fragment check, and determines if the ligands are heavier than a given cutoff
#the user needs to feed in a corresponding file (the smiles string will need to be after the 3rd comma)
#ligands above the cutoff will be printed to a new file with the corresponding line and the MW appended to the end of the line

import os, sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

#get the csv file with smiles
input_file = sys.argv[1]

#get the cutoff mw
cutoff_mw = float(sys.argv[2])

#create the input stream and out stream
#out stream will just append above_#mw to the input name
in_file = open(input_file,"r")

in_file_prefix = input_file.split(".")[0]
out_file_name = in_file_prefix + "_above_" + str(cutoff_mw) + ".csv"
out_file = open(out_file_name,"w")

#read the input file line by line and process
for line in in_file.readlines():
	#get the smiles string and strip off any tailing newline
	smiles_string = line.split(",")[3].strip()
	mol = Chem.MolFromSmiles(smiles_string)

	#get the mass
	mass = 0

	if mol:
		mass = rdMolDescriptors.CalcExactMolWt(mol)

	print(mass,line.strip())

	#if the mass is above the cutoff keep the ligand
	if mass >= cutoff_mw:
		out_file.write(line.strip() + "," + str(mass) + "\n")