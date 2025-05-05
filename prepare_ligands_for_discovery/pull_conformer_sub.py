#helper job for pull_conformers_from_conformer_library.py
#actually pulls down the ligands, extracts the params file, and sets up a test_params directory with all of the ligands for a discovery operation

import os,sys

#arg 1: output location
output_location = sys.argv[1]

#if output location does not end in a "/", make it ned with a "/"
if output_location.endswith("/") == False:
	output_location = output_location + "/"

#arg 2: extract params script
extract_script = sys.argv[2]
#arg 3: fix params file
fix_script = sys.argv[3]

#create the test_params folder
os.system("mkdir test_params")

#create the stream for writing residue types
residue_types_file = open("test_params/residue_types.txt","w")

#write the initial lines in the file
residue_types_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
residue_types_file.write("TYPE_SET_MODE full_atom\n")
residue_types_file.write("ATOM_TYPE_SET fa_standard\n")
residue_types_file.write("ELEMENT_SET default\n")
residue_types_file.write("MM_ATOM_TYPE_SET fa_standard\n")
residue_types_file.write("ORBITAL_TYPE_SET fa_standard\n")
residue_types_file.write("## Test params files\n")
#residue_types_file.write("\n")
#ligand conformer params will be written in following lines

#touch files that are needed to be present, but not necessarily have contents
os.system("touch test_params/patches.txt test_params/exclude_pdb_component_list.txt")

#open ligs.csv file
lig_file = open("ligs.csv", "r")

#hold the previous ligname to avoid redundant downloads and save time
#input list should be sorted with ligane names in alphabetical order, so conformers of the same ligand should be consecutive
#previous_ligname = ""

#hold previous chunk and subchunk instead
previous_chunk = ""
previous_subchunk = ""

#read through each line of the ligs.csv file to extract the conformer param and put in in test_params
for line in lig_file.readlines():
	#line format is: ligname_confnum,chunk,subchunk; extract each component
	ligname = str(line.split(",")[0].split("_")[0])
	confnum = str(line.split(",")[0].split("_")[1])
	chunk = str(line.split(",")[1])
	subchunk = str(line.split(",")[2].strip())

	#print(chunk,previous_chunk)
	#print(subchunk,previous_subchunk)
	#print(chunk != previous_chunk and subchunk != previous_subchunk)

	#delete condensed data from previous run if the current ligname is the same as the previous and previous ligname != ""

	#if ligname != previous_ligname:
	if chunk != previous_chunk or subchunk != previous_subchunk:
		#delete the condensed data to keep workspace clean
		os.system("rm -drf condensed*")
		os.system("s3cmd get  s3://ariosg/ligand_library/" + chunk + "/for_s3/condensed_params_and_db_" + subchunk + ".tar.gz --force --no-progress")

	#set previous ligname to current ligname for next iteration
	#previous_ligname = ligname
	previous_chunk = chunk
	previous_subchunk = subchunk

	#pull the subchunk data from LTS
	#os.system("s3cmd get  s3://ariosg/ligand_library/" + chunk + "/for_s3/condensed_params_and_db_" + subchunk + ".tar.gz")

	#unzip the tar file
	#condensed_params_and_db_0/single_conf_params/Z2354650510_shorthand_params.txt
	os.system("tar -xzf condensed_params_and_db_" + subchunk + ".tar.gz condensed_params_and_db_" + subchunk + "/single_conf_params/" + ligname + "_shorthand_params.txt")

	#go into the folder
	os.chdir("condensed_params_and_db_" + subchunk + "/single_conf_params")

	#run the extraction script
	#python extract_script.py list_file_name.txt confnum ligname
	os.system("python " + extract_script + " " + ligname + "_shorthand_params.txt " + confnum + " " + ligname)

	#run the fix script on the ligand
	os.system("python " + fix_script + " " + ligname + ".params")

	#new file is named with fixed, rename the file
	#also add confnum to the ligand name to make sure that there is no duplicate
	os.system("mv fixed_" + ligname + ".params  " + ligname + "_" + confnum + ".params")

	#move the file to test_params
	os.system("mv " + ligname + "_" + confnum + ".params ../../test_params")

	#write the file name into the residue types file
	residue_types_file.write(ligname + "_" + confnum + ".params\n")

	#move back up
	os.chdir("../..")

#zip the test params
os.system("tar -czf test_params.tar.gz test_params")

#push the test_param directory to the LTS output location
os.system("s3cmd put -r --no-progress test_params" + output_location  +  "test_params.tar.gz")
