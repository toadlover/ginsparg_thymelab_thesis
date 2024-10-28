#The purpose of this script is to look at a directory of placement pdbs and look at its comment data to compose a csv list of placement scores for seeing how Rosetta deals with ligands discovered in other publications versus a random ligand set
#the sores that we want to pull are the system ddg, total motif-like interactions count, total amount of interactions that resemble motifs from the real library, and the ratio of interactions that look real
#this script takes no arguments and is intended to be run in the directory with the placement pdb files


#import
import os,sys

#create a csv file to write to
#in theory, multiple files could be appended together, so do not write a header to the file
#data order to write is: ligand name,placement file name,ddg,interaction count,real motif count, real motif ratio
write_file = open("placement_score_data.csv","w")

#run through the current location and get the data from all pdb files
location = os.getcwd()

for r,d,f in os.walk(location):
	for file in f:
		#only look at pdbs
		if file.endswith(".pdb") and r == location:
			#variables to hold the placement data (besides filename, because we have that)
			ligand_name = ""
			ddg = 0
			interaction_count = 0
			real_motif_count = 0
			real_motif_ratio = 0

			#open file read stream
			file_stream = open(file,"r")

			for line in file_stream.readlines():
				#ligand name
				if line.startswith("Placement: Ligand name:"):
					ligand_name = line.strip().split()[len(line.strip().split()) - 1]

				#ddg
				if line.startswith("Scoring: Post-HighResDock system ddG:"):
					ddg = str(float(line.strip().split()[len(line.strip().split()) - 1]))

				#interaction count
				if line.startswith("Placement motifs: Total motifs made:"):
					interaction_count = str(float(line.strip().split()[len(line.strip().split()) - 1]))				

				#real motif count
				if line.startswith("Placement motifs: Real motif count:"):
					real_motif_count = str(float(line.strip().split()[len(line.strip().split()) - 1]))

				#real motif ratio
				if line.startswith("Placement motifs: Real motif ratio:"):
					real_motif_ratio = str(float(line.strip().split()[len(line.strip().split()) - 1]))

			#once the file is read, make sure we got data and then if we did, write a line to teh csv
			if ligand_name == "":
				continue

			#write the line
			write_file.write(ligand_name + "," + file + "," + ddg + "," + interaction_count + "," + real_motif_count + "," + real_motif_ratio + "\n")

			#test print of line as well:
			print(ligand_name + "," + file + "," + ddg + "," + interaction_count + "," + real_motif_count + "," + real_motif_ratio)
