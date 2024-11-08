#the purpose of this script is to pull placement data csv files from the data bucket and determine where the the paper ligands compare, relative to the library set
#this script utilized s3cmd and assumes that the dat ais stored in a bucket

#import packages
import os,sys

#take in an argument, which is the bucket location to where all relevant csv files for the paper are
#i.e. s3://ariosg/benchmarking_other_papers/ai_powered/ for the AI-powered_Virtual_Screening_of_Large_Compound_Libraries_leads_to_the_Discovery_of_Novel_Inhibitors_of_Sirtuin-1 paper
bucket_path = sys.argv[1]

#make sure that the bucket path ends with a / if it doesn't already
if bucket_path.endswith("/") == False:
	bucket_path = bucket_path + "/"

#first part, collect all csv data for the paper and compile into a new csv file named all_paper_data.csv
#get location of all files via s3cmd ls
#trim the statement to get only the csv files
#there will be 1 file per line
os.system("s3cmd ls -r " + bucket_path + " | grep csv | awk '{print $4}' > csv_file_locations.txt")

#iterate over the csv_file_locations.txt file, pull down each file, and append into the all_paper_data.csv file
locations_file = open("csv_file_locations.txt", "r")

#open a file stream to write the all_paper_data.csv
all_csv = open("all_paper_data.csv", "w")

#write a header line to the csv file
all_csv.write("is_paper_ligand,ligand,ligand_conf,placement_pdb,ddg,interactions,real_motifs_count,real_motifs_ratio")

for line in locations_file.readlines():
	#strip the newline and then get the file
	os.system("s3cmd get " + line.strip())

	#write the file contents to the compiled csv
	#os.system("cat placement_score_data.csv >> all_paper_data.csv")
	#open file stream to read the csv file (we want to add additional data to indicate which are the paper ligands, which is derivable in the path that the file came from)
	small_file = open("placement_score_data.csv", "r")

	#determine whether this is a paper ligand file
	is_paper_ligand = "0"

	if "/paper_ligands/" in line:
		is_paper_ligand = "1"

	#read through the small file and add lines to the all paper file
	#append the paper ligand status to the front of the line
	for data_line in small_file.readlines():

		#pull the first entry of the data_line to extract just the ligand name (everything before the first underscore) so that we can find the best placement for a metric across all conformers
		lig_name = data_line.split("_")[0]

		all_csv.write(is_paper_ligand + "," + lig_name + "," + data_line)

	#delete the file
	os.system("rm placement_score_data.csv")

#close the stream
locations_file.close()
all_csv.close()

#now, read through the compiled csv file and find the best placements for each ligand for ddg and real motifs and put them in a dictionary
#make dictionaries to hold the best of each ligand for ddg, interactions collected, real motif count, and real motif ratio
#the ligand name is used as a key
best_ddg = {}
best_interactions = {}
best_real_count = {}
best_real_ratio = {}

#open the filled out all_paper_data file and read through all the data to potentially add to each dictionary
all_csv = open("all_paper_data.csv","r")

for line in all_csv.readlines():
	#skip the first line in the file
	if line.startswith("is_paper_ligand,ligand,"):
		continue

	#read the line and extract the data and store into a tuple
	#strip the line first to remove the newline
	stripped_line = line.strip()

	#use split to split the line by the data values
	split_line = line.split(",")

	#indices of values in the split:
	#0 - is paper ligand
	#1 - ligand name
	#2 - ligand with conformer
	#3 - placement pdb name
	#4 - ddg
	#5 - interaction count
	#6 - real motif count
	#7 - real motif ratio

	#now go through and check each dictionary to see if the ligand is present and/or if the placement improves upon a previous placement for this ligand

	#ddg
	#if not in dictionary, add as first instance
	if split_line[1] not in best_ddg.keys():
		best_ddg[split_line[1]] = split_line
	#if instance in dictionary, see if new instance is better and replace in dictionary if better
	else:
		if float(best_ddg[split_line[1]][4]) > flost(split_line[4]):
			#replace with better
			best_ddg[split_line[1]] = split_line

	#interaction count
	#if not in dictionary, add as first instance
	if split_line[1] not in best_interactions.keys():
		best_interactions[split_line[1]] = split_line
	#if instance in dictionary, see if new instance is better and replace in dictionary if better
	else:
		if float(best_interactions[split_line[1]][5]) < float(split_line[5]):
			#replace with better
			best_interactions[split_line[1]] = split_line

	#real motif count
	#if not in dictionary, add as first instance
	if split_line[1] not in best_real_count.keys():
		best_real_count[split_line[1]] = split_line
	#if instance in dictionary, see if new instance is better and replace in dictionary if better
	else:
		if float(best_real_count[split_line[1]][6]) < float(split_line[6]):
			#replace with better
			best_real_count[split_line[1]] = split_line

	#real motif ratio
	#if not in dictionary, add as first instance
	if split_line[1] not in best_real_ratio.keys():
		best_real_ratio[split_line[1]] = split_line
	#if instance in dictionary, see if new instance is better and replace in dictionary if better
	else:
		if float(best_real_ratio[split_line[1]][7]) < float(split_line[7]):
			#replace with better
			best_real_ratio[split_line[1]] = split_line

#convert each dictionary to a list so that they can be sorted by relevant values

#ddg
#declare list
best_ddg_list = []
#add each best ddg per ligand to list
for key in best_ddg.keys():
	best_ddg_list.append(best_ddg[key])
#sort list by ddg
best_ddg_list_sorted = sorted(best_ddg_list, key=lambda x: float(x[4]))

#total interactions
#declare list
best_interactions_list = []
#add each best ddg per ligand to list
for key in best_interactions.keys():
	best_interactions_list.append(best_interactions[key])
#sort list by ddg
best_interactions_list_sorted = sorted(best_interactions_list, key=lambda x: float(x[5]), reverse=True)

#real motif count
#declare list
best_real_count_list = []
#add each best ddg per ligand to list
for key in best_real_count.keys():
	best_real_count_list.append(best_real_count[key])
#sort list by ddg
best_real_count_list_sorted = sorted(best_real_count_list, key=lambda x: float(x[6]), reverse=True)

#real motif ratio
#declare list
best_real_ratio_list = []
#add each best ddg per ligand to list
for key in best_real_ratio.keys():
	best_real_ratio_list.append(best_real_ratio[key])
#sort list by ddg
best_real_ratio_list_sorted = sorted(best_real_ratio_list, key=lambda x: float(x[7]), reverse=True)

#now, output each list to csv for data access

#ddg
best_ddg_file = open("best_ddg.csv", "w")
for placement in best_ddg_list_sorted:
	for item in placement:
		best_ddg_file.write(str(item) + ",")
	best_ddg_file.write("\n")

#interaction
best_interactions_file = open("best_interactions.csv", "w")
for placement in best_interactions_list_sorted:
	for item in placement:
		best_interactions_file.write(str(item) + ",")
	best_interactions_file.write("\n")

#real count
best_real_count_file = open("best_real_count.csv", "w")
for placement in best_real_count_list_sorted:
	for item in placement:
		best_real_count_file.write(str(item) + ",")
	best_real_count_file.write("\n")

#real ratio
best_real_ratio_file = open("best_real_ratio.csv", "w")
for placement in best_real_ratio_list_sorted:
	for item in placement:
		best_real_ratio_file.write(str(item) + ",")
	best_real_ratio_file.write("\n")

#finally, run through each list and identify where the paper ligands rank, and output that data
#make a dictionary for the paper ligands to hold the absolute and percentile rank
paper_ligands_ranked = {}

#derive the number of ligands (should only have to measure on one list, as they should all be the same length)
#used to determine the percentile rank of the ligand
num_ligands = len(best_ddg_list_sorted)

#iterate by index so that we also have a counter
#ddg
for i in range(len(best_ddg_list_sorted)):
	#process if the ligand iterated upon is a paper ligand
	if best_ddg_list_sorted[i][0] == "1":
		#add the information for the ligand to the paper_ligands_ranked dictionary
		if best_ddg_list_sorted[i][1] not in paper_ligands_ranked.keys():
			paper_ligands_ranked[best_ddg_list_sorted[i][1]] = [best_ddg_list_sorted[i][1], i + 1, (i + 1) / num_ligands]
		else:
			paper_ligands_ranked[best_ddg_list_sorted[i][1]].append(i + 1)
			paper_ligands_ranked[best_ddg_list_sorted[i][1]].append((i + 1) / num_ligands)

#interactions
for i in range(len(best_interactions_list_sorted)):
	#process if the ligand iterated upon is a paper ligand
	if best_interactions_list_sorted[i][0] == "1":
		#add the information for the ligand to the paper_ligands_ranked dictionary
		if best_interactions_list_sorted[i][1] not in paper_ligands_ranked.keys():
			paper_ligands_ranked[best_interactions_list_sorted[i][1]] = [best_interactions_list_sorted[i][1], i + 1, (i + 1) / num_ligands]
		else:
			paper_ligands_ranked[best_interactions_list_sorted[i][1]].append(i + 1)
			paper_ligands_ranked[best_interactions_list_sorted[i][1]].append((i + 1) / num_ligands)

#real count
for i in range(len(best_real_count_list_sorted)):
	#process if the ligand iterated upon is a paper ligand
	if best_real_count_list_sorted[i][0] == "1":
		#add the information for the ligand to the paper_ligands_ranked dictionary
		if best_real_count_list_sorted[i][1] not in paper_ligands_ranked.keys():
			paper_ligands_ranked[best_real_count_list_sorted[i][1]] = [best_real_count_list_sorted[i][1], i + 1, (i + 1) / num_ligands]
		else:
			paper_ligands_ranked[best_real_count_list_sorted[i][1]].append(i + 1)
			paper_ligands_ranked[best_real_count_list_sorted[i][1]].append((i + 1) / num_ligands)

#real ratio
for i in range(len(best_real_ratio_list_sorted)):
	#process if the ligand iterated upon is a paper ligand
	if best_real_ratio_list_sorted[i][0] == "1":
		#add the information for the ligand to the paper_ligands_ranked dictionary
		if best_real_ratio_list_sorted[i][1] not in paper_ligands_ranked.keys():
			paper_ligands_ranked[best_real_ratio_list_sorted[i][1]] = [best_real_ratio_list_sorted[i][1], i + 1, (i + 1) / num_ligands]
		else:
			paper_ligands_ranked[best_real_ratio_list_sorted[i][1]].append(i + 1)
			paper_ligands_ranked[best_real_ratio_list_sorted[i][1]].append((i + 1) / num_ligands)

#finally, run through paper_ligands_ranked and print out where the paper ligands rank
#make a file to print the ranks
paper_ligands_file = open("paper_ligands_rank.csv", "w")

#write a header for the file
paper_ligands_file.write("ligand,ddg_rank,ddg_percentile,interaction_rank,interaction_percentile,real_motif_count_rank,real_motif_count_percentile,real_motif_ratio_rank,real_motif_ratio_percentile,\n")

for ligand in paper_ligands_ranked.keys():
	for item in ligand:
		paper_ligands_file.write(str(item) + ",")

	paper_ligands_file.write("\n")