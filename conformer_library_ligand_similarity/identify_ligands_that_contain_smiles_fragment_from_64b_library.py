#the purpose of this script is for it to be run on a single short file from the 64b library; there are 
#this script returns a csv list of ligands that contain the input fragment smiles string (agnostic of chirality)
#this looks for .bz files sinze the 64b ligand library is too large to reasonably work with uncompressed

#script arguments:
#1 subchunk library data location
#i.e. /data/project/thymelab/smiles_similarity_of_hits_analysis_space/drug_27/


#3 fragment smiles string
#i.e. CC(C)c1cc(NC(=O)N2Cc3cccc(C#N)c3C2)n(-c2ncccn2)n1

#4 output location for the result list file
#output file name will take the format of ligands_containing_fragment_(subchunk number).csv

#imports
import os,sys
from rdkit import Chem
from rdkit.Chem import AdjustQueryParameters, AdjustQueryProperties, AdjustQueryWhichFlags
import bz2


#read in and process the arguments
library_location = sys.argv[1]

#derive the file name and use it as an extension
extension = library_location.split("/")[len(library_location.split("/")) - 1].split(".bz2")[0]

fragment_smiles = sys.argv[2]


output_location = sys.argv[3]

#modify location to make sure it ends with a /
if output_location.endswith("/") == False:
	output_location = output_location + "/"

#start the output file
output_file = open(output_location + "ligands_containing_fragment_" + str(extension) + ".csv", "w")

#load in the fragment as smiles
fragment_molecule = Chem.MolFromSmiles(fragment_smiles)

"""
#extract query parameters
params = AdjustQueryParameters()
#set makedummiesqueries to true to help make matching more flexible, as matches will be lost otherwise
params.makeDummiesQueries = True
params.aromatizeIfPossible = True

#set this parameter to the fragment molecule
fragment_molecule = AdjustQueryProperties(fragment_molecule, params)

"""

#remove chirality
Chem.RemoveStereochemistry(fragment_molecule)

#make SMARTS of the fragment molecule to pass into the substruct match
fragment_molecule_smarts = Chem.MolToSmarts(fragment_molecule)

fragment_query = Chem.MolFromSmarts(fragment_molecule_smarts)



#open the library file
read_file = open(library_location, "r")

#iterate over the compressed bz2 file
with bz2.open(library_location, 'rt') as file:
	for line in file:
		#extract the smiles string, which is at index 0
		lig = str(line.strip().split()[0])

		#convert the string to smiles
		lig_smiles = Chem.MolFromSmiles(lig)

		#sanity check to make sure that the ligand is ok
		if lig_smiles is None:
			print("Ligand SMILES " + lig + " is invalid")
			continue

		#remove chirality
		Chem.RemoveStereochemistry(lig_smiles)

		#add the makedummies params to the full ligand being investigated
		#lig_smiles = AdjustQueryProperties(lig_smiles, params)

		#run the match of the fragment in the full ligand
		if lig_smiles.HasSubstructMatch(fragment_query, useChirality=False, useQueryQueryMatches=True):
			#if true, the fragment exists within

			#get the ligand name code from the database
			#the name code is either the split index 1 or 2, depending on if stereochemistry is listed
			lig_code = str(line.strip().split()[1])
			if "|" in lig_code:
				lig_code = str(line.strip().split()[2])

			#write the ligand to the output file
			output_file.write(str(line.strip().split()[0]) + "," + lig_code + "," + library_location + "\n")

			print(str(line.strip().split()[0]) + "," + lig_code + "," + library_location + "\n")
