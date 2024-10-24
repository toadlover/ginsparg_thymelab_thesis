#the goal of this script is to create arg files and submit rosetta jobs for rosetta's ligand_discovery_search protocol
#the args file needs to be premade except for the -params_directory_path line, which will be generated by this program
#to run this program, the user needs the working location, the args template file, and the path to the rosetta build
#This script uses slurm to run all discovery jobs in parallel

import os,sys,getpass

#getthe username for use with slurm queue accession
username = getpass.getuser()

#get the directory to work in
working_location = sys.argv[1]

#prepared args file name (and path if needed)
args_file = sys.argv[2]

#rosetta ligand_discovery_Search executable with path
#within rosetta, executable is found at /rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease
#i.e. /scratch/abgvg9/rosetta_may_for_checkin/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease
rosetta_executable = sys.argv[3]

#move to working location
os.chdir(working_location)

#run through each directory in the working location
for r,d,f in os.walk(working_location):
	for dire in d:
		if r == working_location:
			
			#move into the directory
			os.chdir(dire)

			#make a copy of the args file
			os.system("cp " + args_file + " " + dire + "_args")

			#append the path to test params to the args file
			os.system("echo \"-params_directory_path " + r + "/" + dire + "/test_params/\" >> " + dire + "_args")

			#make a slurm job file to run rosetta and then submit
			#make medium partition for 3000 minutes to be safe
			#create and run job_file.slurm
			job_file = open("discover.job", "w")
			job_file.write("#!/bin/bash\n")
			job_file.write("#SBATCH -p medium # Partition to submit to\n")
			job_file.write("#SBATCH -n 1 # Number of cores requested\n")
			job_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
			job_file.write("#SBATCH -t 3000 # Runtime in minutes\n")
			job_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
			job_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
			job_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
			job_file.write(rosetta_executable + " @" + dire + "_args" + "\n")
			job_file.close()

			#sleep for one second and then submit the job
			os.system("sleep 1")
			os.system("sbatch discover.job")

			#check the slurm queue and make sure that the slurm queue is not over 400 jobs (to avoid flooding the queue)
			os.system("squeue -A " + username + " | wc -l  > squeue_file.txt")
			squeue_length_file = open("squeue_file.txt", 'r')
			squeue_length = ""
			for line in squeue_length_file.readlines():
				#rebuild squeue_length to ensure it is only numbers
				for char in line:
					if char.isnumeric():
						squeue_length = squeue_length + char
			squeue_length = int(squeue_length)

			squeue_length_file.close()

			#continually probe for squeue length, and keep cycling until queue dips below 400
			while(squeue_length > 400):
				os.system("squeue -A " + username + " | wc -l  > squeue_file.txt")
				squeue_length_file = open("squeue_file.txt", 'r')
				squeue_length = ""
				for line in squeue_length_file.readlines():
					#rebuild squeue_length to ensure it is only numbers
					for char in line:
						if char.isnumeric():
							squeue_length = squeue_length + char
				squeue_length = int(squeue_length)
				squeue_length_file.close()

			#clean up the slurm queue length file
			os.system("squeue_file.txt")
			#move back up at end
			os.chdir("..")