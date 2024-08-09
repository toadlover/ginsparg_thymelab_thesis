#the purpose of this script is to look at a directory containing folders of conformers in sdf format (1 conformer per file) and create test_params directories for rosetta discovery
import os,sys

#get the directory to work in
working_location = sys.argv[1]

#get the location of the molfile_to_params.py script in the Rosetta package
#it is located at /rosetta/source/scripts/python/public/molfile_to_params.py within Rosetta
#call the entire filename, i.e. /scratch/abgvg9/rosetta_may_for_checkin/rosetta/source/scripts/python/public/molfile_to_params.py
molfile_to_params_path = sys.argv[2]

#move to working location
os.chdir(working_location)

#run through each directory in the working location
for r,d,f in os.walk(working_location):
	for dire in d:
		if r == working_location:
			
			#move into the directory
			os.chdir(dire)

			#make the test_params directory and fill it with params file and set up other nexessary files
			os.system("rm -drf test_params")
			os.system("mkdir test_params")

			#set up residue_types.txt file
			res_types_file = open("test_params/residue_types.txt", "w")
			res_types_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
			res_types_file.write("TYPE_SET_MODE full_atom\n")
			res_types_file.write("ATOM_TYPE_SET fa_standard\n")
			res_types_file.write("ELEMENT_SET default\n")
			res_types_file.write("MM_ATOM_TYPE_SET fa_standard\n")
			res_types_file.write("ORBITAL_TYPE_SET fa_standard\n")
			res_types_file.write("## Test params files\n")
			#res_types_file.write("\n")

			#iterate over each sdf file in the directory and make a params
			for r2,d2,f2 in os.walk(os.getcwd()):
				for file in f2:
					if file.endswith(".sdf"):
						#make the params file
						os.system("python " + molfile_to_params_path + " " + file + " -n " + file.split(".")[0] + " --keep-names --long-names --clobber --no-pdb")

						#write the file to the residue_types file
						res_types_file.write(file + "\n")

						#move the file to test_params
						os.system("mv " + file.split(".")[0] + ".params test_params")

			#touch additional text files
			os.system("touch test_params/exclude_pdb_component_list.txt test_params/patches.txt")

			#close the write stream
			res_types_file.close()

			os.chdir("..")