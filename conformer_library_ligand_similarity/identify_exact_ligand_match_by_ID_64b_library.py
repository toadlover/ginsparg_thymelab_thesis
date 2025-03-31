#far more simple script that just searches the ligand library for a specific ligand name

#script arguments:
#1 subchunk library data location
#i.e. /data/project/thymelab/smiles_similarity_of_hits_analysis_space/drug_27/


#3 full smiles string
#i.e. CC(C)c1cc(NC(=O)N2Cc3cccc(C#N)c3C2)n(-c2ncccn2)n1

#4 output location for the result list file
#output file name will take the format of ligands_containing_fragment_(subchunk number).csv

#imports
import os,sys

import bz2


#read in and process the arguments
library_location = sys.argv[1]

#derive the file name and use it as an extension
extension = library_location.split("/")[len(library_location.split("/")) - 1].split(".bz2")[0]

ligand_id = sys.argv[2]


output_location = sys.argv[3]

#modify location to make sure it ends with a /
if output_location.endswith("/") == False:
	output_location = output_location + "/"

#start the output file
output_file = open(output_location + "ligands_containing_" + str(extension) + ".csv", "w")

#open the library file
read_file = open(library_location, "r")

#iterate over the compressed bz2 file
with bz2.open(library_location, 'rt') as file:
	for line in file:
		
		#if the ligand id is in the line, write it to the out file
		if ligand_id in line:
			#get the ligand name code from the database
			#the name code is either the split index 1 or 2, depending on if stereochemistry is listed
			lig_code = str(line.strip().split()[1])
			if "|" in lig_code:
				lig_code = str(line.strip().split()[2])

			#write the ligand to the output file
			output_file.write(str(line.strip().split()[0]) +  "," + lig_code + "," + library_location + "\n")

			print(str(line.strip().split()[0]) + "," + lig_code + "," + library_location + "\n")
