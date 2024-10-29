#The purpose of this script is to look at directories in the benchmarking_other_papers section of the repository and prepare a file for condor to use to locate where chunks of ligands are to run rosetta on as well as identify where to submit the results
#this looks for files named: condensed_params_and_db_0_output_locations.txt  condensed_params_and_db_1_output_locations.txt  paper_ligands_location.txt
#This script is meant to be run in a benchmarking directory, such as: /ginsparg_thymelab_thesis/benchmarking_other_papers/papers/AI-powered_Virtual_Screening_of_Large_Compound_Libraries_leads_to_the_Discovery_of_Novel_Inhibitors_of_Sirtuin-1

#import
import os,sys

#get current location
location = os.getcwd()

#open a write stream to open a new file to write the file to be read in by the condor job
write_file = open("file_location_data.csv", "w")

#run through the current location and look for the needed files and read them
for r,d,f in os.walk(location):
	for file in f:
		#make sure the file is local
		if r != location:
			continue

		#process if it is one of the 3 desired file names
		if file == "condensed_params_and_db_0_output_locations.txt" or file == "condensed_params_and_db_1_output_locations.txt" or file == "paper_ligands_location.txt":
			print("On file: " + file)

			#open the file stream
			read_stream = open(file, "r")

			#contents of a line should look something like this example:
			#s3://ariosg/benchmarking_other_papers/ai_powered/ligand_inputs/condensed_params_and_db_0/0/test_params.tar.gz
			#we want to extract everything between s3://ariosg/benchmarking_other_papers/ and /test_params.tar.gz; i.e. ai_powered/ligand_inputs/condensed_params_and_db_0/0
			#we want to convert the example string to this: ai_powered/ligand_inputs/condensed_params_and_db_0/0,ai_powered/output_data/condensed_params_and_db_0/0

			#read through the read stream
			for line in read_stream.readlines():
				#strip the newline off the line
				line_no_newline = line.strip()

				#extract the main data we want to keep
				main_string = line_no_newline.split("s3://ariosg/benchmarking_other_papers/")[1].split("/test_params.tar.gz")[0]

				#we now have the first string we want for the input location; now to make the string for output
				#get paper location, which is before the first /
				paper_locator = main_string.split("/")[0]
				#get the unique extension, which is after /ligand_inputs/, which we want to replace with /output_data/
				unique_extension = main_string.split("/ligand_inputs/")[1]

				out_string = paper_locator + "/output_data/" + unique_extension

				#write the main string and out string to the write file
				write_file.write(main_string + "," + out_string + "\n")

			#close the read stream
			read_stream.close()
#close the write stream
write_file.close()

