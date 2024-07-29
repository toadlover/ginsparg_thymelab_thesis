import os,sys

#first argument is the scores list csv that has the placed ligand files with a path leading to them (can include a path)
scores_file = sys.argv[1]

#second argument is the destination folder name (can have a path leading up to it)
destination_folder = sys.argv[2]

#read the scores file
read_file = open(scores_file, "r")

for line in read_file.readlines():
	#break up the line by commas, the file is in positions 0
	placement_file = line.split(",")[0]

	#if the data is "file" (first line), skip
	if placement_file == "file":
		continue

	#otherwise, copy the file to the destination
	os.system("cp " + placement_file + " " + destination_folder)