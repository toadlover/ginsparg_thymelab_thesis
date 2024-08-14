#the purpose of this script is to look at the working location, and push all directories that have a test_params directory to a given location in a bucket
#this program uses the s3cmd command and a slurm job scheduler
#this also uses slurm

import os,sys

#get the directory to work in
working_location = sys.argv[1]

#bucket address in s3 format
#i.e. s3://ariosg/Ox_truncations_12M_chunk_sorted/
bucket_address = sys.argv[2]

#make sure that the bucket address ends with a backslash
if bucket_address.endswith("/") == False:
	bucket_address = bucket_address + "/"

location = working_location

os.chdir(location)

for r,d,f in os.walk(location):
	for dire in d:
		#only looking at immediate directories
		if r == location:
			print(dire)

			#make sure that the test_params directory exists within folder
			has_test_params = False
			for r2,d2,f2 in os.walk(dire):
				for dire2 in d2:
					if dire2 == "test_params":
						has_test_params = True

			#continue of there isn't a test_params directory in the working directory
			if has_test_params == False:
				continue

			#here is where we need to do something.
			#go into the directory, make a job to compress and push, then submit the job
			os.chdir(dire)

			#write a job to compress and put
			write_file = open("compress.job", "w")
			write_file.write("#!/bin/bash\n")
			write_file.write("#SBATCH -p express # Partition to submit to\n")
			write_file.write("#SBATCH -n 1 # Number of cores requested\n")
			write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
			write_file.write("#SBATCH -t 120 # Runtime in minutes\n")
			write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
			write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
			write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
			write_file.write("\n")

			#compress test_params
			write_file.write("tar -czf test_params.tar.gz test_params \n")

			#put the compressed folder in lts
			write_file.write("s3cmd put test_params.tar.gz " + bucket_address + dire + "/ --no-progress \n")

			write_file.close()

			os.system("sleep 0.5")

			#buffer to make sure that squeue does not get overwhelmed
			os.system("squeue -A $(whoami) | wc -l > squeue_length.txt")
			squeue_file = open("squeue_length.txt", "r")
			wait = False
			for line2 in squeue_file.readlines():
				squeue_length = int(line2.strip())
				
				#over 400 lines in file, we need to wait for jobs to finish
				if squeue_length > 450:
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
					if squeue_length < 450:
						wait = False
					else:
						#sleep for 1 second to reduce pressure on slurm
						os.system("sleep 1")
				squeue_file.close()

			#once we overcome the buffer, submit the job and then we can move on
			os.system("sbatch compress.job")

			#move back up at end
			os.chdir("..")
