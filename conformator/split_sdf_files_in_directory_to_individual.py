import os,sys

#get the working location
#get the directory to runconformator on
working_location = sys.argv[1]

#move to working location
os.chdir(working_location)

#run through each sdf file in the working location
for r,d,f in os.walk(working_location):
	for file in f:
		if file.endswith(".sdf") and r == working_location:
			print(file)

			#make a folder based on the name of the ligand (up to the first underscore)
			file_prefix = file.split("_")[0]
			os.system("mkdir " + file_prefix)

			#run obabel to split the ligands
			os.system("obabel " + file + " -O " + file_prefix + "/" + file_prefix + "_.sdf -m")

			#run through each newly made file and adjust the ligand name in the file
			for r2,d2,f2 in os.walk(file_prefix):
				for single_file in f2:
					#read the file and write to a temporary copy
					read_file = open(r2 + "/" + single_file,"r")
					write_file = open(r2 + "/temp.sdf", "w")

					#line counter, we are only interested in line 1
					line_counter = 0

					for line in read_file.readlines():
						if line_counter == 0:
							write_file.write(single_file.split(".")[0] + "\n")
						else:
							write_file.write(line)

						line_counter = line_counter + 1

					read_file.close()
					write_file.close()

					#write the temp over the original
					os.system("mv " r2 + "/temp.sdf " + r2 + "/" + single_file)