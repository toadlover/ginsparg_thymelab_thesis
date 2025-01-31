#the purpose of this script is to look at a directory of directories of params and/or test_params files and send them to a bucket location


#example call: python compress_and_send_test_params_to_bucket.py /data/project/thymelab/smiles_similarity_of_hits_analysis_space/drug_27/drug_27_best_sorted_test_chiral_splits/receiving/0 s3://ariosg/drug27_closest_ligands/

#imports
import os,sys

#arguments
#working location; the head of the directories
#i.e. /data/project/thymelab/smiles_similarity_of_hits_analysis_space/drug_27/drug_27_best_sorted_test_chiral_splits/receiving/0 where the directories immediately below are named after ligands, and each ligand directory contains either compressed test_params or a list of params files
working_location = sys.argv[1]

#location in bucket to send each individual directory to
bucket_location = sys.argv[2]
#if the bucket location does not end with a /, make it end with a /
if bucket_location.endswith("/") == False:
	bucket_location = bucket_location + "/"

#move into the working location
os.chdir(working_location)

#iterate over each directory and push compressed test_params to the directory
for r,d,f in os.walk(working_location):
	for dire in d:
		#enter the directory
		os.chdir(dire)

		#check if there is a prepared test_params directory
		if os.path.exists("test_params.tar.gz"):
			#send the directory up to the bucket
			os.system("s3cmd put test_params.tar.gz " + bucket_location + dire + "/test_params.tar.gz")

		#otherwise, take the existing params, prepare a test_params directory, compress it, and then push it
		else:
			os.system("rm *sdf")
			os.system("mkdir test_params")
			os.system("mv *.params test_params")

			#move into the test_params folder and set it up
			os.chdir("test_params")

			#touch necessary files
			os.system("touch exclude_pdb_component_list.txt patches.txt")

			#make the residue_types file
			res_types_file = open("residue_types.txt", "w")
			res_types_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
			res_types_file.write("TYPE_SET_MODE full_atom\n")
			res_types_file.write("ATOM_TYPE_SET fa_standard\n")
			res_types_file.write("ELEMENT_SET default\n")
			res_types_file.write("MM_ATOM_TYPE_SET fa_standard\n")
			res_types_file.write("ORBITAL_TYPE_SET fa_standard\n")
			res_types_file.write("## Test params files\n")
			#write all params file names to this file
			for r,d,f in os.walk(os.getcwd()):
				for file in f:
					if file.endswith(".params"):
						res_types_file.write(file + "\n")

			res_types_file.close()

			#go back up
			os.chdir("..")

			#compress the directory and toss the original
			os.system("tar -czf  test_params.tar.gz test_params")

			#os.system("rm -drf test_params")
			#send the directory up to the bucket
			os.system("s3cmd put test_params.tar.gz " + bucket_location + dire)

			os.chdir("..")



