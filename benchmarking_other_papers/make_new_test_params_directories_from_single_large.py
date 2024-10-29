#the purpose of this script is to clean up oversized test_params directories (>100 params) and split it into new test_params directories with fewer params (up to 50)
#the script takes in a path of a test_params directory, and then at the level of the directory (at, not in), it makes numbered folders with smaller test params directories inside (with compressed directories as well)

#import 
import os,sys

test_params_location = sys.argv[1]

#make sure that we actually have a test_params directory
if test_params_location.endswith("test_params") == False or est_params_location.endswith("test_params/") == False:
	print("Input of: " + test_params_location + " is bad.")
	quit()

#move to the test_params location and then go one level up
os.chdir(test_params_location)
os.chdir("..")

#make the first new test_params directory
test_params_directory_counter = 0
conf_counter = 0

test_params_directory_counter = test_params_directory_counter + 1
os.system("mkdir " + str(test_params_directory_counter))
os.system("mkdir " + str(test_params_directory_counter) + "/test_params")

#touch files that are needed to be present, but not necessarily have contents
os.system("touch " + str(test_params_directory_counter) + "/test_params/patches.txt " + str(test_params_directory_counter) + "/test_params/exclude_pdb_component_list.txt")

#start making the residue_types_file
residue_types_file = open(str(test_params_directory_counter) + "/test_params/residue_types.txt","w")

#write the initial lines in the file
residue_types_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
residue_types_file.write("TYPE_SET_MODE full_atom\n")
residue_types_file.write("ATOM_TYPE_SET fa_standard\n")
residue_types_file.write("ELEMENT_SET default\n")
residue_types_file.write("MM_ATOM_TYPE_SET fa_standard\n")
residue_types_file.write("ORBITAL_TYPE_SET fa_standard\n")
residue_types_file.write("## Test params files\n")

#look through all files in the original test_params directory and copy them out
for r,d,f in os.walk(test_params_location):
	for file in f:
		#work with params files in the test_params location
		if file.endswith(".params") and r == test_params_location:
			#increment the conformer counter
			conf_counter = conf_counter + 1

			#copy the conformer into the new test_params and write the conformer to the residue_types
			os.system("cp " + r + "/" + file + " " + str(test_params_directory_counter) + "/test_params")
			residue_types_file.write(file + "\n")

			#if we hit 50 conformers, compress the current test_params and then set up the next test_params directory
			if conf_counter % 50 == 0:
				os.chdir(str(test_params_directory_counter))
				os.system("tar -czf test_params.tar.gz test_params")
				os.chdir("..")

				residue_types_file.close()

				test_params_directory_counter = test_params_directory_counter + 1
				os.system("mkdir " + str(test_params_directory_counter))
				os.system("mkdir " + str(test_params_directory_counter) + "/test_params")

				#touch files that are needed to be present, but not necessarily have contents
				os.system("touch " + str(test_params_directory_counter) + "/test_params/patches.txt " + str(test_params_directory_counter) + "/test_params/exclude_pdb_component_list.txt")

				#start making the residue_types_file
				residue_types_file = open(str(test_params_directory_counter) + "/test_params/residue_types.txt","w")

				#write the initial lines in the file
				residue_types_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
				residue_types_file.write("TYPE_SET_MODE full_atom\n")
				residue_types_file.write("ATOM_TYPE_SET fa_standard\n")
				residue_types_file.write("ELEMENT_SET default\n")
				residue_types_file.write("MM_ATOM_TYPE_SET fa_standard\n")
				residue_types_file.write("ORBITAL_TYPE_SET fa_standard\n")
				residue_types_file.write("## Test params files\n")

#finish with final directory
os.chdir(str(test_params_directory_counter))
os.system("tar -czf test_params.tar.gz test_params")
os.chdir("..")

residue_types_file.close()