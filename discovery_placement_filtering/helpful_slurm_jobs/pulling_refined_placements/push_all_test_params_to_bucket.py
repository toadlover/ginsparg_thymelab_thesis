#this is a fairly specific script that may not be useful for other people to use
#the purpose of this script is to look at a given folder, and look down at all test_params.tar.gz folders
#for each folder, push it to a user-specified bucket location. when pushing to the bucket, an extra folder extension will be used that utilizes the ligand name that was holdingthe tar.gz in order to avoid collision and make selection easier
#this also makes a list file that osg will use for queueing up each job based on the ligands

#imports
import os,sys

#get the bucket location
bucket_location = sys.argv[1]

#have bucket location end with a /
if bucket_location.endswith("/") == False:
	bucket_location = bucket_location + "/"

#get the folder location, will look down for ALL test_params from here
working_location = sys.argv[2]

if working_location.endswith("/") == False:
	working_location = working_location + "/"

#open write stream to write ligand names to text file for condor queue text file
write_file = open(working_location + "ligands_list.txt", "w")

for r,d,f in os.walk(working_location):
	for file in f:
		if file == "test_params.tar.gz":
			#attempt to push test params to bucket

			#derive the ligand name as the lowest directory in the file root
			ligand = r.split("/")[len(r.split("/")) - 1]

			write_file.write(ligand + "\n")

			print("s3cmd put " + r + "/test_params.tar.gz " + bucket_location + ligand + "/")