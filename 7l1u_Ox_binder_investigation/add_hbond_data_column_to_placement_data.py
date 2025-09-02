#the purpose of this script is to add a column to an existing placements data file for a list of any hydrogen bond energies
#two columns will be added, one for a colon-separated list of residues with hbond scores with real interactions, and one for non-real interactions

#imports
import os,sys

#read in the input file
#this file assumes that the first column is a pdb placement file with either a path to it (safe), or the file is local (dangerous)
input_file = sys.argv[1]

#make a string for the output file, add _hbonds before the csv extension
output_file = input_file.split(".")[0] + "_hbonds.csv"

read_file = open(input_file,"r")

#create a write stream
write_file = open(output_file,"w")

#read the file
for line in read_file.readlines():
	#create a stripped version of the line
	stripped_line = line.strip()

	#handle the header
	if line.startswith("file,ddg"):
		write_file.write(stripped_line + ",real_hbonds,nonreal_hbonds\n")

		continue

	#otherwise, read the file
	placements_file = line.split(",")[0]

	placements_stream = open(placements_file,"r")

	#create strings for holding real and not real hbonds
	real_hbonds_string = ""
	non_real_hbonds_string = ""

	#read the file
	for line2 in placements_stream.readlines():
		#if "Real motif check" is in the line, we can check the line
		if "Real motif check" in line2:
			#obtain the hbond score, if it is under -0.01, call it an hbond and get the residue associated
			hbond_score = float(line2.split("Hbond_score:")[2].split()[0])
			if hbond_score > -0.01:
				continue

			#derive the residue
			residue = line2.split("Hbond_score:")[1].split("_")[1]

			#combine the residue and score, tail it with a semicolon
			combined_term = residue + ":" + str(hbond_score) + ";"

			#determine if it is a real motif or not
			if "No real match" in line2:
				non_real_hbonds_string = non_real_hbonds_string + combined_term
			else:
				real_hbonds_string = real_hbonds_string + combined_term

	#done with all motifs, add the new data to the columns for the row
	write_file.write(stripped_line + "," + real_hbonds_string + "," + non_real_hbonds_string + "\n")
