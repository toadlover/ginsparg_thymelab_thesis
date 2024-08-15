import os,sys

#untarred sdf file containing up to 5k ligands
starting_file = sys.argv[1]

#license key used to activate conformator
license_key = sys.argv[2]

#split name to remove extension
#if there is a path, remove that too to work in the current location
starting_file_no_ext = starting_file.split("/")[len(starting_file.split("/")) - 1].split(".")[0]

#make a directory to perform all operations in and then move into it
os.system("mkdir " + starting_file_no_ext)
os.system("cp " + starting_file + " " + starting_file_no_ext)
os.chdir(starting_file_no_ext)

#activate conformator
os.system("/conformator_for_container/conformator_1.2.1/conformator --license " + license_key)

#run conformator on the starting file
os.system("/conformator_for_container/conformator_1.2.1/conformator -i " + starting_file_no_ext + ".sdf -o confs.sdf --keep3d --hydrogens -n 15 -v 0")

#use babel to break the confs.sdf file into individual sdf files; also apply unique names to the ligand name (per conformer)
os.system("babel confs.sdf individual_conf.sdf -m")

#make a molecule name dictionary to hold the molecule name and number of times that it has been encountered
molecule_name_dict = {}

#make directories to hold single ligand conformer sdfs and params
os.system("mkdir single_conf_sdfs single_conf_params")

#run through each individual_conf#.sdf, get the molecule name and assign it a unique number identifier, rename the file with the unique name and append the file to a confs file with the unique name
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		if file.startswith("individual_conf") and file.endswith(".sdf"):
			#get the first line of the file, which has the molecule name
			#if this line is blank, continue (we failed somewhere, there should be a name)
			read_file = open(file, "r")

			molecule_name = ""
			for line in read_file.readlines():

				molecule_name = str(line.strip("\n"))

				#add molecule name to molecule name dict if not present (with count of 1), otherwise increment 1 to the count
				if molecule_name not in molecule_name_dict:
					molecule_name_dict[molecule_name] = 1
				else:
					molecule_name_dict[molecule_name] = molecule_name_dict[molecule_name] + 1

				break

			read_file.close()

			#mol name + unique id number
			unique_name = molecule_name + "_" + str(molecule_name_dict[molecule_name])

			#now we have the molecule name and a count for it; re-write the file with the unique number identifier after it in the name (also name the file with it)
			write_file = open(unique_name + ".sdf", "w")

			#open the same read file again so we can copy it
			read_file = open(file, "r")

			#bool to indicate if we are on the first line or not
			first_line = True

			for line in read_file.readlines():
				#first line, write the unique name here
				if first_line == True:
					first_line = False
					write_file.write(unique_name + "\n")
				#copy the line from the read file
				else:
					write_file.write(line)

			read_file.close()
			write_file.close()

			#append the written file with the updated name to a master confs list
			#this will be used for later use with shapedb/VAMS
			os.system("cat " + unique_name + ".sdf >> confs_named.sdf")

			#append the unique name to a file list for being able to access file names when we do stuff with shape.db/VAMS
			os.system("echo " + unique_name + " >> " + starting_file_no_ext + "_lig_name_list.txt")

			#make a params file of the unique file
			os.system("python /conformator_for_container/molfile_to_params.py " + unique_name + ".sdf -n " + unique_name + " --keep-names --long-names --clobber --no-pdb")

			#at this point, we should have made a sdf and params file for the single molecule file
			#move the sdf and param to their own folders
			os.system("mv " + unique_name + ".sdf single_conf_sdfs")
			os.system("mv " + unique_name + ".params single_conf_params")

			#delete the individual_conf#.sdf since we don't need it and have a better copy with the unique name sdf
			os.system("rm " + file)

#all sdf and params files should be made
#remove confs.sdf since we have build confs_named.sdf
os.system("rm confs.sdf")

os.chdir("..")

#tar confs_named, single_conf_sdfs, and single_conf_params
#move back up a directory and tar starting_file_no_ext; delete the original
os.system("tar -czvf " + starting_file_no_ext + ".tar.gz " + starting_file_no_ext)
os.system("rm -dr " + starting_file_no_ext)


"""
os.system("tar -czvf confs_named.sdf.tar.gz confs_named.sdf")
os.system("rm confs_named.sdf")
os.system("tar -czvf single_conf_sdfs.tar.gz single_conf_sdfs")
os.system("rm -dr single_conf_sdfs")
os.system("tar -czvf single_conf_params.tar.gz single_conf_params")
os.system("rm -dr single_conf_params")
"""
