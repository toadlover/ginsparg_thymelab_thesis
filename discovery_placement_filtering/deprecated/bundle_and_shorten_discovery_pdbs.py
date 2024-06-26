import os,sys

#counter to make directories to break things up into bundles of 100 files
counter = 0

#make directory
os.system("mkdir " + str(counter))

#open a file to write the file names in
write_file = open(str(counter) + "/original_file_names.txt", "w")

#sub counter to count the individual ligands
sub_counter = 0
print(counter)
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		if file.endswith(".pdb"):
			#extract the ligand name
			lig_name = file.split("_Trio")[1].split("_motif")[0].split("_")[1]
			conf_num = file.split("_Trio")[1].split("_motif")[0].split("_")[2]
			lig_conf = lig_name + "_" + conf_num

			#copy the current file as ligconf with the sub counter into the counter directory
			os.system("cp " + file + " " + str(counter) + "/" + lig_conf + "_" + str(sub_counter) + ".pdb")

			#write the original file name to the write file
			write_file.write(file + " " + str(counter) + "/" + lig_conf + "_" + str(sub_counter) + ".pdb\n")

			#increment the counters and behavior if hitting 100
			sub_counter = sub_counter + 1

			if sub_counter % 100 == 0:
				
				counter = counter + 1
				print(counter)
				os.system("mkdir " + str(counter))
				write_file.close()
				#open a file to write the file names in
				write_file = open(str(counter) + "/original_file_names.txt", "w")
