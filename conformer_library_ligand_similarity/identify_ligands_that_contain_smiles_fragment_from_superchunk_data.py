#the purpose of this script is for it to be run on a single superchunk of the ligand library
#This script is to be run on pulled down data that is in csv form that was a result from run_smiles_comparison_check_on_superchunk_of_library.py, with the ligand data (including a smiles string)
#this script returns a csv list of ligands that contain the input fragment smiles string (agnostic of chirality)

#example file line from input superchunk csv file:
#40740,7,Z2475112899,0.12244897959183673,CC(C)c1cc(NC(=O)N2Cc3cccc(C#N)c3C2)n(-c2ncccn2)n1

#script arguments:
#1 subchunk library data location
#i.e. /data/project/thymelab/smiles_similarity_of_hits_analysis_space/drug_27/

#2 subchunk number
#value between 0-531

#3 fragment smiles string
#i.e. CC(C)c1cc(NC(=O)N2Cc3cccc(C#N)c3C2)n(-c2ncccn2)n1

#4 output location for the result list file
#output file name will take the format of ligands_containing_fragment_(subchunk number).csv

#imports
import os,sys
from rdkit import Chem
from rdkit.Chem import AdjustQueryParameters, AdjustQueryProperties, AdjustQueryWhichFlags


#read in and process the arguments
library_location = sys.argv[1]

#modify location to make sure it ends with a /
if library_location.endswith("/") == False:
	library_location = library_location + "/"

subchunk = sys.argv[2]
fragment_smiles = sys.argv[3]


output_location = sys.argv[4]

#modify location to make sure it ends with a /
if output_location.endswith("/") == False:
	output_location = output_location + "/"

#start the output file
output_file = open(output_location + "ligands_containing_fragment_" + str(subchunk) + ".csv", "w")

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

#make SMARTS of the fragment molecule to pass into the substruct match
fragment_molecule_smarts = Chem.MolToSmarts(fragment_molecule)

fragment_query = Chem.MolFromSmarts(fragment_molecule_smarts)

#remove chirality
Chem.RemoveStereochemistry(fragment_molecule)

#iterate over all files in the subchunk
for r,d,f in os.walk(library_location + str(subchunk)):
	for file in f:
		#only look at relevant csv files
		if file.endswith("_smiles_similarity.csv"):
			
			print(r + "/" + file)

			#read the file and see if the fragment exists in each ligand smiles
			read_file = open(r + "/" + file, "r")

			for line in read_file.readlines():
				#extract the smiles string
				lig = str(line.strip().split(",")[4])

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
					#write the ligand to the output file
					output_file.write(str(line.strip().split(",")[0]) + "," + str(line.strip().split(",")[1]) + "," + str(line.strip().split(",")[2]) + "," + str(line.strip().split(",")[4]) + "\n")
					print(str(line.strip().split(",")[0]) + "," + str(line.strip().split(",")[1]) + "," + str(line.strip().split(",")[2]) + "," + str(line.strip().split(",")[4]))
