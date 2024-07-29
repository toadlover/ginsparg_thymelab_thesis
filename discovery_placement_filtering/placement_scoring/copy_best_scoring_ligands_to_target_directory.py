import os,sys

# python /data/user/abgvg9/ginsparg_thymelab_thesis/discovery_placement_filtering/placement_scoring/copy_best_scoring_ligands_to_target_directory.py  /scratch/abgvg9/discovery_results/top_1000_placement/12M_raw_ddg_sorted.csv /scratch/abgvg9/discovery_results/top_1000_placement/agonist_12M_passing_placements/

#first argument is the scores list csv that has the placed ligand files with a path leading to them (can include a path)
scores_file = sys.argv[1]

#second argument is the destination folder name (can have a path leading up to it)
destination_folder = sys.argv[2]

#append a backslash if there is not one currently
if destination_folder.endswith("/") == False:
	destination_folder = destination_folder + "/"

#read the scores file
read_file = open(scores_file, "r")

#break up the files into sections of up to 100 files by separating into directories so as to not overwhelm pymol if these files are converted into sessions
sub_dir = 0
#make directory
os.system("mkdir " + destination_folder + str(sub_dir))

file_counter = 0

for line in read_file.readlines():
	#break up the line by commas, the file is in positions 0
	placement_file = line.split(",")[0]

	#if the data is "file" (first line), skip
	if placement_file == "file":
		continue

	#increment the file counter
	file_counter = file_counter + 1

	#if the file counter is divisible by 100, incremend the sub_dir
	if file_counter % 100 == 0:
		sub_dir = sub_dir + 1
		#make directory
		os.system("mkdir " + destination_folder + str(sub_dir))

	#otherwise, copy the file to the destination
	os.system("cp " + placement_file + " " + destination_folder + str(sub_dir))