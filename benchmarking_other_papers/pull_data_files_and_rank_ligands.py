#the purpose of this script is to pull placement data csv files from the data bucket and determine where the the paper ligands compare, relative to the library set
#this script utilized s3cmd and assumes that the dat ais stored in a bucket

#import packages
import os,sys

#take in an argument, which is the bucket location to where all relevant csv files for the paper are
#i.e. s3://ariosg/benchmarking_other_papers/ai_powered/ for the AI-powered_Virtual_Screening_of_Large_Compound_Libraries_leads_to_the_Discovery_of_Novel_Inhibitors_of_Sirtuin-1 paper
bucket_path = sys.argv[1]

#first part, collect all csv data for the paper and compile into a new csv file named all_paper_data.csv
#get location of all files via s3cmd ls
#trim the statement to get only the csv files
#there will be 1 file per line
os.system("s3cmd ls -r s3://ariosg/benchmarking_other_papers/ | grep csv | awk '{print $4}' > csv_file_locations.txt")

#iterate over the csv_file_locations.txt file, pull down each file, and append into the all_paper_data.csv file
locations_file = open("csv_file_locations.txt", "r")

#open a file stream to write the all_paper_data.csv
all_csv = open("all_paper_data.csv", "w")

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
		all_csv.write(is_paper_ligand + "," + data_line)

	#delete the file
	os.system("rm placement_score_data.csv")

#close the stream
locations_file.close()
all_csv.close()

#now, read through the compiled csv file and find the best placements for each ligand for ddg and real motifs and list them
#make lists