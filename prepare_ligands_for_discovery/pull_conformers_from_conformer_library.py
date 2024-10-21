import os,sys

#the sourse csv file with all conformers
#listed in the format:
#ligname_confnum,chunk,subchunk
source_csv = sys.argv[1]

#the exact path to output chunks of ligands to
#follow with a slash
#i.e. s3://ariosg/Ox_truncations_12M/
s3_output_dir = sys.argv[2]

if s3_output_dir.endswith("/") == False:
	s3_output_dir = s3_output_dir + "/"

#the number of ligands per chunk subsection of output library
#trying 1000 for now, but want this to be changeable
#will make a job to set up each chunk of the x ligands
ligands_per_section = int(sys.argv[3])

#location to the sub script to be called upon by the jobs
#pull_conformer_sub.py
location_of_sub = sys.argv[4]

#paths to the scripts to extract single ligand and fix the file
#ginsparg_thymelab_thesis/params_file_compression/extract_single_param_from_condensed_file.py
#ginsparg_thymelab_thesis/params_file_compression/fix_condensed_param_file_spacing.py
extract_file = sys.argv[5]
fix_file = sys.argv[6]

#open the source csv
source_csv_file = open(source_csv,"r")

#line counter
line_counter = 0

#file counter
file_counter = 0

#create a main directory to hold output operations and a subdirectory for the first file of output operations
os.system("mkdir operations")

os.system("mkdir operations/0")

#open a file to start writing lines to in the operations section
small_csv = open("operations/0/ligs.csv", "w")

#read the source csv and prepare to fire off jobs for every ligands per section
for line in source_csv_file.readlines():
	
	#increment the counter
	line_counter = line_counter + 1

	#write the line to the small csv
	small_csv.write(line)

	#if the line counter % ligands per section = 0, submit a job to process the small csv to set it up for discovery
	if line_counter % ligands_per_section == 0:
		
		#close the small csv write stream
		small_csv.close()

		#write a job to process the ligands in the file
		write_file = open("operations/" + str(file_counter) + "/process.job","w")
		write_file.write("#!/bin/bash\n")
		write_file.write("#SBATCH -p short # Partition to submit to\n")
		write_file.write("#SBATCH -n 1 # Number of cores requested\n")
		write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
		write_file.write("#SBATCH -t 720 # Runtime in minutes\n")
		write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
		write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
		write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
		write_file.write("\n")
		#python sub_file.py  s3_output_dir/file_counter #this is the output location for these conformers to be pulled for discovery
		write_file.write("python " + location_of_sub + " " + s3_output_dir + str(file_counter) + " " + extract_file + " " + fix_file + " \n")
		write_file.close()


		#submit the job
		os.chdir("operations/" + str(file_counter))
		os.system("sbatch process.job")
		os.chdir("../..")

		#sleep 0.5 seconds (to ensure no overload of slurm queue)
		os.system("sleep 0.5")

		#make sure that we do not overwhelm the slurm queue, have personal queue be no longer than 400 lines
		#pull slurm queue length
		os.system("squeue -A $(whoami) | wc -l > squeue_length.txt")
		squeue_file = open("squeue_length.txt", "r")
		wait = False
		for line2 in squeue_file.readlines():
			squeue_length = int(line2.strip())
			
			#over 400 lines in file, we need to wait for jobs to finish
			if squeue_length > 400:
				wait = True
				#sleep for 1 second to reduce pressure on slurm
				os.system("sleep 1")
		squeue_file.close()
		while wait:
			#probe the queue again and keep probing until squeue_length is < 400
			os.system("squeue -A $(whoami) | wc -l > squeue_length.txt")
			squeue_file = open("squeue_length.txt", "r")
			for line2 in squeue_file.readlines():
				squeue_length = int(line2.strip())
				
				#over 400 lines in file, we need to wait for jobs to finish
				if squeue_length < 400:
					wait = False
				else:
					#sleep for 1 second to reduce pressure on slurm
					os.system("sleep 1")
			squeue_file.close()




		#increment the file counter and make a new folder for the next round of outputs
		file_counter = file_counter + 1

		os.system("mkdir operations/" + str(file_counter))

		#make new output write stream
		small_csv = open("operations/" + str(file_counter) + "/ligs.csv", "w")

#final submission of remaining ligs.csv


#write a job to process the ligands in the file
write_file = open("operations/" + str(file_counter) + "/process.job","w")
write_file.write("#!/bin/bash\n")
write_file.write("#SBATCH -p express # Partition to submit to\n")
write_file.write("#SBATCH -n 1 # Number of cores requested\n")
write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
write_file.write("#SBATCH -t 120 # Runtime in minutes\n")
write_file.write("#SBATCH --mem=3000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
write_file.write("\n")
#python sub_file.py  s3_output_dir/file_counter 
#this is the output location for these conformers to be pulled for discovery
write_file.write("python " + location_of_sub + " " + s3_output_dir + str(file_counter) + "  \n")
write_file.close()


#submit the job
os.chdir("operations/" + str(file_counter))
os.system("sbatch process.job")
os.chdir("../..")