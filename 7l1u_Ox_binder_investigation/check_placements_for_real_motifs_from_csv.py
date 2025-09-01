#import os,sys

#the purpose is to read through a raw scores csvfile, and determine if placements in the file have read motifs against residues of interest
#for this specific case, we want to have real motif interactions with ASN20 and at least one of ASP211,GLU217,ARG328,HIS335 (using indexing found in the motifs section)

#the user will give the csv file in the script call so it can be read over, and an outputted filtered csv will be returned

import os,sys

csv_file = sys.argv[1]

#create read and write streams
read_csv = open(csv_file, "r")

write_csv = open("filtered_" + csv_file, "w")

#have a line counter to help track how many files have been processed
line_counter = 0

#read the csv file
for line in read_csv.readlines():
	line_counter = line_counter + 1
	#print progress every 10,000 files processed
	if line_counter % 10000 == 0:
		print("Processed: " + str(line_counter))

	#handle the header
	if line.startswith("file,ddg,"):
		write_csv.write(line)
		continue

	#otherwise, handle a file (which should have its path included, for easy reading)
	placement_file = line.split(",")[0]
	read_file = open(placement_file,"r")

	#declare variables to hold whether the placement file has interactions with the peptide and receptor residues of interest, default as false and change true
	has_peptide_interaction = False
	has_receptor_interaction = False

	for line2 in read_file.readlines():
		#check if the line is a real motifs check line
		if line2.startswith("Placement motifs: Real motif check"):
			#if so, check and see if it has a match (if it has "No real match", then pass)
			if "No real match" in line2:
				continue

			#check and see if any of the residues of interest are involved
			if "_ASN20_" in line2:
				has_peptide_interaction = True

			#ASP211,GLU217,ARG328,HIS335
			if "_ASP211_" in line2 or "_GLU217_" in line2 or "_ARG328_" in line2 or "_HIS335_" in line2:
				has_receptor_interaction = True

	#we have now processed the file, determine whether to print hte line or not (print the line if both values are true)
	if has_receptor_interaction and has_peptide_interaction:
		write_csv.write(line)



