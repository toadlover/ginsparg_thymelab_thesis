#The purpose of this script is to work in the current location and prepare an entire compressed subchunk of ligands into a set of test_params directories for use with Rosetta discovery
#This script is made in mind with being used for the benchmarking Rosetta on the discovery of other papers, and used to make a set of ligands with similar MW to ligands from the paper
#This script will put all conformers of all ligands in the subchunk into a test_params directories and will make directories of up to 100 conformers
#This script takes an input subchunk location in the bucket (and will pull down with s3cmd), and an output location for all discovery directories

#example call: python /data/user/abgvg9/ginsparg_thymelab_thesis/prepare_ligands_for_discovery/prepare_whole_subchunk_for_discovery.py s3://ariosg/ligand_library/00297/for_s3/condensed_params_and_db_0.tar.gz s3://ariosg/benchmarking_other_papers/ai_powered/ligand_inputs/

#initial import
import os,sys

#retain the start of the working location
starting_location = os.getcwd()

#get the subchunk with location and the bucket output location

#i.e. s3://ariosg/ligand_library/00297/for_s3/condensed_params_and_db_0.tar.gz
subchunk_location = sys.argv[1]

#extract the subchunk name
#strip path and tar extension
#s3://ariosg/ligand_library/00297/for_s3/condensed_params_and_db_0.tar.gz -> condensed_params_and_db_0
subchunk_name = subchunk_location.split("/")[len(subchunk_location.split("/")) - 1].split(".tar")[0]

#i.e. s3://ariosg/benchmarking_other_papers/ai_powered/ligand_inputs/
output_location = sys.argv[2]

#if the output location does not end with a /, make it end with a /
#outputs will be compressed test_params directories so that they take up less space and be more easily pulled for OSG condor
if output_location.endswith("/") == False:
	output_location = output_location + "/"

#add the subchunk name to the output location, in case more subchunks are wanted and then we won't overwrite a previous different chunk
output_location = output_location + subchunk_name + "/"

#pull the subchunk to the working location
os.system("s3cmd get " + subchunk_location)

#unzip the subchunk
os.system("tar -xzf " + subchunk_name + ".tar.gz")

#derive paths to the params extraction script and params fix script
extract_script = os.path.dirname(os.path.abspath(__file__)) + "/../params_file_compression/extract_single_param_from_condensed_file.py"
fix_script = os.path.dirname(os.path.abspath(__file__)) + "/../params_file_compression/fix_condensed_param_file_spacing.py"

#count the number of conformers processed
#use to get a final count and to know when to make the next test_params_directory
conformer_counter = 0
test_params_directory_count = 0

#make a stream to write a file that has all of the output test_params directories for tracking and could be used for later jobs
output_list = open(subchunk_name + "_output_locations.txt","w")

#start making a test params directory
os.system("mkdir test_params")

#touch files that are needed to be present, but not necessarily have contents
os.system("touch test_params/patches.txt test_params/exclude_pdb_component_list.txt")

#start making the residue_types_file
residue_types_file = open("test_params/residue_types.txt","w")

#write the initial lines in the file
residue_types_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
residue_types_file.write("TYPE_SET_MODE full_atom\n")
residue_types_file.write("ATOM_TYPE_SET fa_standard\n")
residue_types_file.write("ELEMENT_SET default\n")
residue_types_file.write("MM_ATOM_TYPE_SET fa_standard\n")
residue_types_file.write("ORBITAL_TYPE_SET fa_standard\n")
residue_types_file.write("## Test params files\n")





#run through each shorthand params file under chunkname/single_conf_params
for r,d,f in os.walk(subchunk_name + "/single_conf_params"):
	#run through each file and work with the shorthand_params files
	for file in f:
		#make sure it is a shorthand file
		if file.endswith("_shorthand_params.txt"):
			#read the file and determine how many conformers there are to derive
			shorthand_file = open(r + "/" + file,"r")

			#hold the conformer number
			conf_num = 0

			for line in shorthand_file.readlines():
				#we want the line that starts with _a:, since we can derive the number of conformers for this ligand from it
				if line.startswith("_a:"):
					conf_num = len(line.split(",")) - 1
					break

			#extract the ligand name for naming params files
			lig_name = file.split("_")[0]

			#now, iterate over each ligand and extract it into the test_params
			#range from 1 - #confs
			for i in range(1,conf_num + 1):
				#string for the ligand name with conformer (no extension)
				ligconf = lig_name + "_" + str(i)

				#run the extract script on the condensed file
				os.system("python " + extract_script + " " + r + "/" + file + " " + str(i) + " " + ligconf)

				#run the fix script on the conformer
				os.system("python " + fix_script + " " + ligconf + ".params")

				#move the fixed file to the test_params directory
				os.system("mv fixed_" + ligconf + ".params test_params")

				#write the conformer to the residue_types file
				residue_types_file.write(ligconf + "\n")

				#increment the conformer counter
				conformer_counter = conformer_counter + 1

				#if we hit 100 ligands (counter % 100 == 0), we are done with this test_params
				#compress it and send to the output bucket and then prepare for a new test_params
				if conformer_counter % 100 == 0:
					#close the residue_types write stream
					residue_types_file.close()

					#compress test_params
					os.system("tar -czf test_params.tar.gz test_params")

					#push the directory to the bucket
					os.system("s3cmd put test_params.tar.gz " + output_location + "/" + str(test_params_directory_count) + "/test_params.tar.gz")

					#write the test_params location to be outputted to
					output_list.write(output_location + str(test_params_directory_count) + "/test_params.tar.gz\n")

					#increment the test_params counter
					test_params_directory_count = test_params_directory_count + 1

					#reset the test_params directory
					os.system("rm test_params/*params test_params/residue_types.txt")

					#start a new residue_types write stream
					residue_types_file = open("test_params/residue_types.txt","w")

					#write the initial lines in the file
					residue_types_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
					residue_types_file.write("TYPE_SET_MODE full_atom\n")
					residue_types_file.write("ATOM_TYPE_SET fa_standard\n")
					residue_types_file.write("ELEMENT_SET default\n")
					residue_types_file.write("MM_ATOM_TYPE_SET fa_standard\n")
					residue_types_file.write("ORBITAL_TYPE_SET fa_standard\n")
					residue_types_file.write("## Test params files\n")


#at the end, wrap up the final test_params directory beign worked on
#close the residue_types write stream
residue_types_file.close()

#compress test_params
os.system("tar -czf test_params.tar.gz test_params")

#push the directory to the bucket
os.system("s3cmd put test_params.tar.gz " + output_location + "/" + str(test_params_directory_count) + "/test_params.tar.gz")

#write the test_params location to be outputted to
output_list.write(output_location + str(test_params_directory_count) + "/test_params.tar.gz")




