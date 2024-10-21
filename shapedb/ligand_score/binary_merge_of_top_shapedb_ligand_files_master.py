#second step in selecting top shape-scoring ligands
#use this after the list of the top hits per super chunk are pulled by get_top_ligands_in_chunks_sub.py
#this will go through each resulting file for the superchunk (since they are generated in parallel), and look at all of the resulting files to pull the best of the best
#this operates and uses get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py to merge 2 files of the starting 531 superchunk top x ligand files by shapedb score for a given shape
#this program controls the use of that python script to pair 2 files together to merge until the data condenses all superchunks together
#the program utilizes slurm to run the binary merges in parallel
#each merge takes the x ligands in each file, combines them together (2x), and then uses a heap and takes the top x of the 2x

import os,sys

#going to try to use a heap data structure to better keep the best scores
import heapq

#get args which are: the file prefix, and max number of ligands to keep

#file prefix i.e. suvo, OxB_7_shifted; do not have trailing underscore unless there is one (full file name is prefix_NN_subchunk_chunk.tar.gz)
file_prefix = str(sys.argv[1])

#value to hold the number of ligands to keep (ideally the same number as the max number in the files being looked at, but allowing the option to change if desired)
#currently using 8,000,000 when looking at 8 truncated shapes of OxA/OxB
max_ligands_to_keep = int(sys.argv[2])

#get the path to where get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py is (in case we are not working in ginsparg_thymelab_thesis)
path_to_top_ligand_script = str(sys.argv[3])

#if the path does not end with a /, make it end with a /
if path_to_top_ligand_script.endswith("/") == False:
	path_to_top_ligand_script = path_to_top_ligand_script + "/"

#get your login for use with slurm queue monitoring
username = os.getlogin()

#log(2)531 is a little larger than 9 (necessetates 9-10 rounds of binary merging)
#there is a smarter iterative way to do this, but for now I will write all 9 steps out and maybe clean it up later

#create a list to carry the file names down the script and from the prior step
#address will overwrite prior in subsequent steps
address_names = []
prior_address_names = []

#create a list to hold all active slurm jobs so that we can check the slurm queue for any active jobs submitted on the current step
active_slurm_jobs = []

#first round of pairing of 531 (x) files; 530/2 = 265, merge 265 times, last file will be carried down to next step
#last file address is 53000_53099
for i in range(0,265):
	#determine the corresponding superchunks to write the string for based on i and then derive an identifier name for the merged file
	f1_min = i * 200
	f1_max = f1_min + 99
	f2_min = f1_min + 100
	f2_max = f2_min + 99

	#make string
	f1_min = str(f1_min)
	while len(f1_min) < 5:
		f1_min = "0" + f1_min
	f1_max = str(f1_max)
	while len(f1_max) < 5:
		f1_max = "0" + f1_max
	f2_min = str(f2_min)
	while len(f2_min) < 5:
		f2_min = "0" + f2_min
	f2_max = str(f2_max)
	while len(f2_max) < 5:
		f2_max = "0" + f2_max
	print(f1_min + "_" + f1_max,f2_min + "_" + f2_max)

	#make addresses for the inputs and the output
	f1_address = f1_min + "_" + f1_max
	f2_address = f2_min + "_" + f2_max
	new_address = f1_min + "_" + f2_max

	address_names.append(new_address)

	#write a slurm job to binary merge 2 files
	write_file = open("binary_merge.job","w")
	write_file.write("#!/bin/bash\n")
	write_file.write("#SBATCH -p express # Partition to submit to\n")
	write_file.write("#SBATCH -n 1 # Number of cores requested\n")
	write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	write_file.write("#SBATCH -t 20 # Runtime in minutes\n")
	write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
	write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	write_file.write("\n")
	#python get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py OxB_7_shifted 8000000 45400_45499 45500_45599 merged
	write_file.write("python " + path_to_top_ligand_script + "get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py " + file_prefix + " " + str(max_ligands_to_keep) + " " + f1_address + " " + f2_address + " " + new_address + "  \n")
	write_file.close()

	#submit the file and collect the slurm job id
	os.system("sbatch binary_merge.job > slurm_job_id.txt")

	#read the slurm job id file and collect the id so it can be added to the slurm id list
	id_file = open("slurm_job_id.txt", "r")
	for idline in id_file.readlines():
		jobid = idline.split()[3].strip()
		active_slurm_jobs.append(jobid)
		break

#add 53000_53099 to end of address names so it gets used in the next round, there should be 266 entries at this point
address_names.append("53000_53099")

#all round 1 jobs should be submitted, check slurm queue until the queue no longer has jobs in the active slurm jobs list
jobs_running = True
while jobs_running:
	#set job running to False, can turn back to true
	jobs_running = False

	#pull the slurm queue 
	#if anyone else ever uses this, this line needs to change (or the smarter thing would be to use whoami; I should probably do that...)
	os.system("squeue -A " + username + " > slurm_queue.txt")

	#buffer in a 1 second sleep to make sure we don't ping slurm too frequently
	os.system("sleep 1")

	#read the slurm queue
	slurm_queue = open("slurm_queue.txt", "r")
	for slurm_line in slurm_queue.readlines():
		#extract the first entry in the line (which would be the id), and determine if it is in the list of active slurm jobs
		#if it is, set jobs running to true, because this id is in the list of active jobs
		cur_id = slurm_line.split()[0]

		if cur_id in active_slurm_jobs:
			jobs_running = True

#----------------
#at this point, we should be through the first run, now condense 266 files into 133

#clear the active slurm jobs list
active_slurm_jobs = []

#set addressnames to prior and empty addressnames
prior_address_names = address_names
address_names = []

#run through everything in the prior address names and do binary merge on each
for i in range(0,int(len(prior_address_names)/2)):
	#determine the corresponding superchunks to write the string for based on i and then derive an identifier name for the merged file

	f1_address = prior_address_names[2 * i]
	f2_address = prior_address_names[(2 * i) + 1]

	#generate the new address by taking the min of f1 address and max of f2 address

	new_address = f1_address.split("_")[0] + "_" + f2_address.split("_")[1]

	address_names.append(new_address)

	#write a slurm job to binary merge 2 files
	write_file = open("binary_merge.job","w")
	write_file.write("#!/bin/bash\n")
	write_file.write("#SBATCH -p express # Partition to submit to\n")
	write_file.write("#SBATCH -n 1 # Number of cores requested\n")
	write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	write_file.write("#SBATCH -t 20 # Runtime in minutes\n")
	write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
	write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	write_file.write("\n")
	#python get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py OxB_7_shifted 8000000 45400_45499 45500_45599 merged
	write_file.write("python " + path_to_top_ligand_script + "get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py " + file_prefix + " " + str(max_ligands_to_keep) + " " + f1_address + " " + f2_address + " " + new_address + "  \n")
	write_file.close()

	#submit the file and collect the slurm job id
	os.system("sbatch binary_merge.job > slurm_job_id.txt")

	#read the slurm job id file and collect the id so it can be added to the slurm id list
	id_file = open("slurm_job_id.txt", "r")
	for idline in id_file.readlines():
		jobid = idline.split()[3].strip()
		active_slurm_jobs.append(jobid)
		break

#check if there is an even or odd number of entries in prior address names. If there is, add the last entry to address names for the next round
#this is easier than thinking through what the name of the skipped odd address is if there even is one for the round
if len(prior_address_names) % 2 == 1:
	address_names.append(prior_address_names[len(prior_address_names) - 1])

#all current round jobs jobs should be submitted, check slurm queue until the queue no longer has jobs in the active slurm jobs list
jobs_running = True
while jobs_running:
	#set job running to False, can turn back to true
	jobs_running = False

	#pull the slurm queue 
	#if anyone else ever uses this, this line needs to change (or the smarter thing would be to use whoami; I should probably do that...)
	os.system("squeue -A " + username + " > slurm_queue.txt")

	#buffer in a 1 second sleep to make sure we don't ping slurm too frequently
	os.system("sleep 1")

	#read the slurm queue
	slurm_queue = open("slurm_queue.txt", "r")
	for slurm_line in slurm_queue.readlines():
		#extract the first entry in the line (which would be the id), and determine if it is in the list of active slurm jobs
		#if it is, set jobs running to true, because this id is in the list of active jobs
		cur_id = slurm_line.split()[0]

		if cur_id in active_slurm_jobs:
			jobs_running = True

#----------------
#at this point, we should be through the second run, now condense 133 files into 66 (with 1 leftover to pass down)

#clear the active slurm jobs list
active_slurm_jobs = []

#set addressnames to prior and empty addressnames
prior_address_names = address_names
address_names = []

#run through everything in the prior address names and do binary merge on each
for i in range(0,int(len(prior_address_names)/2)):
	#determine the corresponding superchunks to write the string for based on i and then derive an identifier name for the merged file

	f1_address = prior_address_names[2 * i]
	f2_address = prior_address_names[(2 * i) + 1]

	#generate the new address by taking the min of f1 address and max of f2 address

	new_address = f1_address.split("_")[0] + "_" + f2_address.split("_")[1]

	address_names.append(new_address)

	#write a slurm job to binary merge 2 files
	write_file = open("binary_merge.job","w")
	write_file.write("#!/bin/bash\n")
	write_file.write("#SBATCH -p express # Partition to submit to\n")
	write_file.write("#SBATCH -n 1 # Number of cores requested\n")
	write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	write_file.write("#SBATCH -t 20 # Runtime in minutes\n")
	write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
	write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	write_file.write("\n")
	#python get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py OxB_7_shifted 8000000 45400_45499 45500_45599 merged
	write_file.write("python " + path_to_top_ligand_script + "get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py " + file_prefix + " " + str(max_ligands_to_keep) + " " + f1_address + " " + f2_address + " " + new_address + "  \n")
	write_file.close()

	#submit the file and collect the slurm job id
	os.system("sbatch binary_merge.job > slurm_job_id.txt")

	#read the slurm job id file and collect the id so it can be added to the slurm id list
	id_file = open("slurm_job_id.txt", "r")
	for idline in id_file.readlines():
		jobid = idline.split()[3].strip()
		active_slurm_jobs.append(jobid)
		break

#check if there is an even or odd number of entries in prior address names. If there is, add the last entry to address names for the next round
#this is easier than thinking through what the name of the skipped odd address is if there even is one for the round
if len(prior_address_names) % 2 == 1:
	address_names.append(prior_address_names[len(prior_address_names) - 1])

#all current round jobs jobs should be submitted, check slurm queue until the queue no longer has jobs in the active slurm jobs list
jobs_running = True
while jobs_running:
	#set job running to False, can turn back to true
	jobs_running = False

	#pull the slurm queue 
	#if anyone else ever uses this, this line needs to change (or the smarter thing would be to use whoami; I should probably do that...)
	os.system("squeue -A " + username + " > slurm_queue.txt")

	#buffer in a 1 second sleep to make sure we don't ping slurm too frequently
	os.system("sleep 1")

	#read the slurm queue
	slurm_queue = open("slurm_queue.txt", "r")
	for slurm_line in slurm_queue.readlines():
		#extract the first entry in the line (which would be the id), and determine if it is in the list of active slurm jobs
		#if it is, set jobs running to true, because this id is in the list of active jobs
		cur_id = slurm_line.split()[0]

		if cur_id in active_slurm_jobs:
			jobs_running = True

#----------------
#at this point, we should be through the third run, now condense 67 files into 33 (with 1 leftover to pass down)

#clear the active slurm jobs list
active_slurm_jobs = []

#set addressnames to prior and empty addressnames
prior_address_names = address_names
address_names = []

#run through everything in the prior address names and do binary merge on each
for i in range(0,int(len(prior_address_names)/2)):
	#determine the corresponding superchunks to write the string for based on i and then derive an identifier name for the merged file

	f1_address = prior_address_names[2 * i]
	f2_address = prior_address_names[(2 * i) + 1]

	#generate the new address by taking the min of f1 address and max of f2 address

	new_address = f1_address.split("_")[0] + "_" + f2_address.split("_")[1]

	address_names.append(new_address)

	#write a slurm job to binary merge 2 files
	write_file = open("binary_merge.job","w")
	write_file.write("#!/bin/bash\n")
	write_file.write("#SBATCH -p express # Partition to submit to\n")
	write_file.write("#SBATCH -n 1 # Number of cores requested\n")
	write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	write_file.write("#SBATCH -t 20 # Runtime in minutes\n")
	write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
	write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	write_file.write("\n")
	#python get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py OxB_7_shifted 8000000 45400_45499 45500_45599 merged
	write_file.write("python " + path_to_top_ligand_script + "get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py " + file_prefix + " " + str(max_ligands_to_keep) + " " + f1_address + " " + f2_address + " " + new_address + "  \n")
	write_file.close()

	#submit the file and collect the slurm job id
	os.system("sbatch binary_merge.job > slurm_job_id.txt")

	#read the slurm job id file and collect the id so it can be added to the slurm id list
	id_file = open("slurm_job_id.txt", "r")
	for idline in id_file.readlines():
		jobid = idline.split()[3].strip()
		active_slurm_jobs.append(jobid)
		break

#check if there is an even or odd number of entries in prior address names. If there is, add the last entry to address names for the next round
#this is easier than thinking through what the name of the skipped odd address is if there even is one for the round
if len(prior_address_names) % 2 == 1:
	address_names.append(prior_address_names[len(prior_address_names) - 1])

#all current round jobs jobs should be submitted, check slurm queue until the queue no longer has jobs in the active slurm jobs list
jobs_running = True
while jobs_running:
	#set job running to False, can turn back to true
	jobs_running = False

	#pull the slurm queue 
	#if anyone else ever uses this, this line needs to change (or the smarter thing would be to use whoami; I should probably do that...)
	os.system("squeue -A " + username + " > slurm_queue.txt")

	#buffer in a 1 second sleep to make sure we don't ping slurm too frequently
	os.system("sleep 1")

	#read the slurm queue
	slurm_queue = open("slurm_queue.txt", "r")
	for slurm_line in slurm_queue.readlines():
		#extract the first entry in the line (which would be the id), and determine if it is in the list of active slurm jobs
		#if it is, set jobs running to true, because this id is in the list of active jobs
		cur_id = slurm_line.split()[0]

		if cur_id in active_slurm_jobs:
			jobs_running = True

#----------------
#at this point, we should be through the fourth run, now condense 34 files into 17 (with 0 leftover to pass down)

#clear the active slurm jobs list
active_slurm_jobs = []

#set addressnames to prior and empty addressnames
prior_address_names = address_names
address_names = []

#run through everything in the prior address names and do binary merge on each
for i in range(0,int(len(prior_address_names)/2)):
	#determine the corresponding superchunks to write the string for based on i and then derive an identifier name for the merged file

	f1_address = prior_address_names[2 * i]
	f2_address = prior_address_names[(2 * i) + 1]

	#generate the new address by taking the min of f1 address and max of f2 address

	new_address = f1_address.split("_")[0] + "_" + f2_address.split("_")[1]

	address_names.append(new_address)

	#write a slurm job to binary merge 2 files
	write_file = open("binary_merge.job","w")
	write_file.write("#!/bin/bash\n")
	write_file.write("#SBATCH -p express # Partition to submit to\n")
	write_file.write("#SBATCH -n 1 # Number of cores requested\n")
	write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	write_file.write("#SBATCH -t 20 # Runtime in minutes\n")
	write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
	write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	write_file.write("\n")
	#python get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py OxB_7_shifted 8000000 45400_45499 45500_45599 merged
	write_file.write("python " + path_to_top_ligand_script + "get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py " + file_prefix + " " + str(max_ligands_to_keep) + " " + f1_address + " " + f2_address + " " + new_address + "  \n")
	write_file.close()

	#submit the file and collect the slurm job id
	os.system("sbatch binary_merge.job > slurm_job_id.txt")

	#read the slurm job id file and collect the id so it can be added to the slurm id list
	id_file = open("slurm_job_id.txt", "r")
	for idline in id_file.readlines():
		jobid = idline.split()[3].strip()
		active_slurm_jobs.append(jobid)
		break

#check if there is an even or odd number of entries in prior address names. If there is, add the last entry to address names for the next round
#this is easier than thinking through what the name of the skipped odd address is if there even is one for the round
if len(prior_address_names) % 2 == 1:
	address_names.append(prior_address_names[len(prior_address_names) - 1])

#all current round jobs jobs should be submitted, check slurm queue until the queue no longer has jobs in the active slurm jobs list
jobs_running = True
while jobs_running:
	#set job running to False, can turn back to true
	jobs_running = False

	#pull the slurm queue 
	#if anyone else ever uses this, this line needs to change (or the smarter thing would be to use whoami; I should probably do that...)
	os.system("squeue -A " + username + " > slurm_queue.txt")

	#buffer in a 1 second sleep to make sure we don't ping slurm too frequently
	os.system("sleep 1")

	#read the slurm queue
	slurm_queue = open("slurm_queue.txt", "r")
	for slurm_line in slurm_queue.readlines():
		#extract the first entry in the line (which would be the id), and determine if it is in the list of active slurm jobs
		#if it is, set jobs running to true, because this id is in the list of active jobs
		cur_id = slurm_line.split()[0]

		if cur_id in active_slurm_jobs:
			jobs_running = True

#----------------
#at this point, we should be through the fifth run, now condense 17 files into 8 (with 1 leftover to pass down)

#clear the active slurm jobs list
active_slurm_jobs = []

#set addressnames to prior and empty addressnames
prior_address_names = address_names
address_names = []

#run through everything in the prior address names and do binary merge on each
for i in range(0,int(len(prior_address_names)/2)):
	#determine the corresponding superchunks to write the string for based on i and then derive an identifier name for the merged file

	f1_address = prior_address_names[2 * i]
	f2_address = prior_address_names[(2 * i) + 1]

	#generate the new address by taking the min of f1 address and max of f2 address

	new_address = f1_address.split("_")[0] + "_" + f2_address.split("_")[1]

	address_names.append(new_address)

	#write a slurm job to binary merge 2 files
	write_file = open("binary_merge.job","w")
	write_file.write("#!/bin/bash\n")
	write_file.write("#SBATCH -p express # Partition to submit to\n")
	write_file.write("#SBATCH -n 1 # Number of cores requested\n")
	write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	write_file.write("#SBATCH -t 20 # Runtime in minutes\n")
	write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
	write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	write_file.write("\n")
	#python get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py OxB_7_shifted 8000000 45400_45499 45500_45599 merged
	write_file.write("python " + path_to_top_ligand_script + "get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py " + file_prefix + " " + str(max_ligands_to_keep) + " " + f1_address + " " + f2_address + " " + new_address + "  \n")
	write_file.close()

	#submit the file and collect the slurm job id
	os.system("sbatch binary_merge.job > slurm_job_id.txt")

	#read the slurm job id file and collect the id so it can be added to the slurm id list
	id_file = open("slurm_job_id.txt", "r")
	for idline in id_file.readlines():
		jobid = idline.split()[3].strip()
		active_slurm_jobs.append(jobid)
		break

#check if there is an even or odd number of entries in prior address names. If there is, add the last entry to address names for the next round
#this is easier than thinking through what the name of the skipped odd address is if there even is one for the round
if len(prior_address_names) % 2 == 1:
	address_names.append(prior_address_names[len(prior_address_names) - 1])

#all current round jobs jobs should be submitted, check slurm queue until the queue no longer has jobs in the active slurm jobs list
jobs_running = True
while jobs_running:
	#set job running to False, can turn back to true
	jobs_running = False

	#pull the slurm queue 
	#if anyone else ever uses this, this line needs to change (or the smarter thing would be to use whoami; I should probably do that...)
	os.system("squeue -A " + username + " > slurm_queue.txt")

	#buffer in a 1 second sleep to make sure we don't ping slurm too frequently
	os.system("sleep 1")

	#read the slurm queue
	slurm_queue = open("slurm_queue.txt", "r")
	for slurm_line in slurm_queue.readlines():
		#extract the first entry in the line (which would be the id), and determine if it is in the list of active slurm jobs
		#if it is, set jobs running to true, because this id is in the list of active jobs
		cur_id = slurm_line.split()[0]

		if cur_id in active_slurm_jobs:
			jobs_running = True

#----------------
#at this point, we should be through the sixth run, now condense 9 files into 4 (with 1 leftover to pass down)

#clear the active slurm jobs list
active_slurm_jobs = []

#set addressnames to prior and empty addressnames
prior_address_names = address_names
address_names = []

#run through everything in the prior address names and do binary merge on each
for i in range(0,int(len(prior_address_names)/2)):
	#determine the corresponding superchunks to write the string for based on i and then derive an identifier name for the merged file

	f1_address = prior_address_names[2 * i]
	f2_address = prior_address_names[(2 * i) + 1]

	#generate the new address by taking the min of f1 address and max of f2 address

	new_address = f1_address.split("_")[0] + "_" + f2_address.split("_")[1]

	address_names.append(new_address)

	#write a slurm job to binary merge 2 files
	write_file = open("binary_merge.job","w")
	write_file.write("#!/bin/bash\n")
	write_file.write("#SBATCH -p express # Partition to submit to\n")
	write_file.write("#SBATCH -n 1 # Number of cores requested\n")
	write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	write_file.write("#SBATCH -t 20 # Runtime in minutes\n")
	write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
	write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	write_file.write("\n")
	#python get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py OxB_7_shifted 8000000 45400_45499 45500_45599 merged
	write_file.write("python " + path_to_top_ligand_script + "get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py " + file_prefix + " " + str(max_ligands_to_keep) + " " + f1_address + " " + f2_address + " " + new_address + "  \n")
	write_file.close()

	#submit the file and collect the slurm job id
	os.system("sbatch binary_merge.job > slurm_job_id.txt")

	#read the slurm job id file and collect the id so it can be added to the slurm id list
	id_file = open("slurm_job_id.txt", "r")
	for idline in id_file.readlines():
		jobid = idline.split()[3].strip()
		active_slurm_jobs.append(jobid)
		break

#check if there is an even or odd number of entries in prior address names. If there is, add the last entry to address names for the next round
#this is easier than thinking through what the name of the skipped odd address is if there even is one for the round
if len(prior_address_names) % 2 == 1:
	address_names.append(prior_address_names[len(prior_address_names) - 1])

#all current round jobs jobs should be submitted, check slurm queue until the queue no longer has jobs in the active slurm jobs list
jobs_running = True
while jobs_running:
	#set job running to False, can turn back to true
	jobs_running = False

	#pull the slurm queue 
	#if anyone else ever uses this, this line needs to change (or the smarter thing would be to use whoami; I should probably do that...)
	os.system("squeue -A " + username + " > slurm_queue.txt")

	#buffer in a 1 second sleep to make sure we don't ping slurm too frequently
	os.system("sleep 1")

	#read the slurm queue
	slurm_queue = open("slurm_queue.txt", "r")
	for slurm_line in slurm_queue.readlines():
		#extract the first entry in the line (which would be the id), and determine if it is in the list of active slurm jobs
		#if it is, set jobs running to true, because this id is in the list of active jobs
		cur_id = slurm_line.split()[0]

		if cur_id in active_slurm_jobs:
			jobs_running = True

#----------------
#at this point, we should be through the seventh run, now condense 5 files into 2 (with 1 leftover to pass down)

#clear the active slurm jobs list
active_slurm_jobs = []

#set addressnames to prior and empty addressnames
prior_address_names = address_names
address_names = []

#run through everything in the prior address names and do binary merge on each
for i in range(0,int(len(prior_address_names)/2)):
	#determine the corresponding superchunks to write the string for based on i and then derive an identifier name for the merged file

	f1_address = prior_address_names[2 * i]
	f2_address = prior_address_names[(2 * i) + 1]

	#generate the new address by taking the min of f1 address and max of f2 address

	new_address = f1_address.split("_")[0] + "_" + f2_address.split("_")[1]

	address_names.append(new_address)

	#write a slurm job to binary merge 2 files
	write_file = open("binary_merge.job","w")
	write_file.write("#!/bin/bash\n")
	write_file.write("#SBATCH -p express # Partition to submit to\n")
	write_file.write("#SBATCH -n 1 # Number of cores requested\n")
	write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	write_file.write("#SBATCH -t 20 # Runtime in minutes\n")
	write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
	write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	write_file.write("\n")
	#python get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py OxB_7_shifted 8000000 45400_45499 45500_45599 merged
	write_file.write("python " + path_to_top_ligand_script + "get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py " + file_prefix + " " + str(max_ligands_to_keep) + " " + f1_address + " " + f2_address + " " + new_address + "  \n")
	write_file.close()

	#submit the file and collect the slurm job id
	os.system("sbatch binary_merge.job > slurm_job_id.txt")

	#read the slurm job id file and collect the id so it can be added to the slurm id list
	id_file = open("slurm_job_id.txt", "r")
	for idline in id_file.readlines():
		jobid = idline.split()[3].strip()
		active_slurm_jobs.append(jobid)
		break

#check if there is an even or odd number of entries in prior address names. If there is, add the last entry to address names for the next round
#this is easier than thinking through what the name of the skipped odd address is if there even is one for the round
if len(prior_address_names) % 2 == 1:
	address_names.append(prior_address_names[len(prior_address_names) - 1])

#all current round jobs jobs should be submitted, check slurm queue until the queue no longer has jobs in the active slurm jobs list
jobs_running = True
while jobs_running:
	#set job running to False, can turn back to true
	jobs_running = False

	#pull the slurm queue 
	#if anyone else ever uses this, this line needs to change (or the smarter thing would be to use whoami; I should probably do that...)
	os.system("squeue -A " + username + " > slurm_queue.txt")

	#buffer in a 1 second sleep to make sure we don't ping slurm too frequently
	os.system("sleep 1")

	#read the slurm queue
	slurm_queue = open("slurm_queue.txt", "r")
	for slurm_line in slurm_queue.readlines():
		#extract the first entry in the line (which would be the id), and determine if it is in the list of active slurm jobs
		#if it is, set jobs running to true, because this id is in the list of active jobs
		cur_id = slurm_line.split()[0]

		if cur_id in active_slurm_jobs:
			jobs_running = True

#----------------
#at this point, we should be through the eighth run, now condense 3 files into 1 (with 1 leftover to pass down)

#clear the active slurm jobs list
active_slurm_jobs = []

#set addressnames to prior and empty addressnames
prior_address_names = address_names
address_names = []

#run through everything in the prior address names and do binary merge on each
for i in range(0,int(len(prior_address_names)/2)):
	#determine the corresponding superchunks to write the string for based on i and then derive an identifier name for the merged file

	f1_address = prior_address_names[2 * i]
	f2_address = prior_address_names[(2 * i) + 1]

	#generate the new address by taking the min of f1 address and max of f2 address

	new_address = f1_address.split("_")[0] + "_" + f2_address.split("_")[1]

	address_names.append(new_address)

	#write a slurm job to binary merge 2 files
	write_file = open("binary_merge.job","w")
	write_file.write("#!/bin/bash\n")
	write_file.write("#SBATCH -p express # Partition to submit to\n")
	write_file.write("#SBATCH -n 1 # Number of cores requested\n")
	write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	write_file.write("#SBATCH -t 20 # Runtime in minutes\n")
	write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
	write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	write_file.write("\n")
	#python get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py OxB_7_shifted 8000000 45400_45499 45500_45599 merged
	write_file.write("python " + path_to_top_ligand_script + "get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py " + file_prefix + " " + str(max_ligands_to_keep) + " " + f1_address + " " + f2_address + " " + new_address + "  \n")
	write_file.close()

	#submit the file and collect the slurm job id
	os.system("sbatch binary_merge.job > slurm_job_id.txt")

	#read the slurm job id file and collect the id so it can be added to the slurm id list
	id_file = open("slurm_job_id.txt", "r")
	for idline in id_file.readlines():
		jobid = idline.split()[3].strip()
		active_slurm_jobs.append(jobid)
		break

#check if there is an even or odd number of entries in prior address names. If there is, add the last entry to address names for the next round
#this is easier than thinking through what the name of the skipped odd address is if there even is one for the round
if len(prior_address_names) % 2 == 1:
	address_names.append(prior_address_names[len(prior_address_names) - 1])

#all current round jobs jobs should be submitted, check slurm queue until the queue no longer has jobs in the active slurm jobs list
jobs_running = True
while jobs_running:
	#set job running to False, can turn back to true
	jobs_running = False

	#pull the slurm queue 
	#if anyone else ever uses this, this line needs to change (or the smarter thing would be to use whoami; I should probably do that...)
	os.system("squeue -A " + username + " > slurm_queue.txt")

	#buffer in a 1 second sleep to make sure we don't ping slurm too frequently
	os.system("sleep 1")

	#read the slurm queue
	slurm_queue = open("slurm_queue.txt", "r")
	for slurm_line in slurm_queue.readlines():
		#extract the first entry in the line (which would be the id), and determine if it is in the list of active slurm jobs
		#if it is, set jobs running to true, because this id is in the list of active jobs
		cur_id = slurm_line.split()[0]

		if cur_id in active_slurm_jobs:
			jobs_running = True

#----------------
#at this point, we should be through the ninth and final run, now condense 2 files into 1 (with 1 leftover to pass down)

#clear the active slurm jobs list
active_slurm_jobs = []

#set addressnames to prior and empty addressnames
prior_address_names = address_names
address_names = []

#run through everything in the prior address names and do binary merge on each
for i in range(0,int(len(prior_address_names)/2)):
	#determine the corresponding superchunks to write the string for based on i and then derive an identifier name for the merged file

	f1_address = prior_address_names[2 * i]
	f2_address = prior_address_names[(2 * i) + 1]

	#generate the new address by taking the min of f1 address and max of f2 address

	new_address = f1_address.split("_")[0] + "_" + f2_address.split("_")[1]

	address_names.append(new_address)

	#write a slurm job to binary merge 2 files
	write_file = open("binary_merge.job","w")
	write_file.write("#!/bin/bash\n")
	write_file.write("#SBATCH -p express # Partition to submit to\n")
	write_file.write("#SBATCH -n 1 # Number of cores requested\n")
	write_file.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	write_file.write("#SBATCH -t 20 # Runtime in minutes\n")
	write_file.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	write_file.write("#SBATCH -o hostname_%A_%a.out # Standard out goes to this file\n")
	write_file.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	write_file.write("\n")
	#python get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py OxB_7_shifted 8000000 45400_45499 45500_45599 merged
	write_file.write("python " + path_to_top_ligand_script + "get_top_ligands_top_chunk_top_lists_between_2_super_chunk_lists.py " + file_prefix + " " + str(max_ligands_to_keep) + " " + f1_address + " " + f2_address + " " + new_address + "  \n")
	write_file.close()

	#submit the file and collect the slurm job id
	os.system("sbatch binary_merge.job > slurm_job_id.txt")

	#read the slurm job id file and collect the id so it can be added to the slurm id list
	id_file = open("slurm_job_id.txt", "r")
	for idline in id_file.readlines():
		jobid = idline.split()[3].strip()
		active_slurm_jobs.append(jobid)
		break

#check if there is an even or odd number of entries in prior address names. If there is, add the last entry to address names for the next round
#this is easier than thinking through what the name of the skipped odd address is if there even is one for the round
if len(prior_address_names) % 2 == 1:
	address_names.append(prior_address_names[len(prior_address_names) - 1])

#all current round jobs jobs should be submitted, check slurm queue until the queue no longer has jobs in the active slurm jobs list
jobs_running = True
while jobs_running:
	#set job running to False, can turn back to true
	jobs_running = False

	#pull the slurm queue 
	#if anyone else ever uses this, this line needs to change (or the smarter thing would be to use whoami; I should probably do that...)
	os.system("squeue -A " + username + " > slurm_queue.txt")

	#buffer in a 1 second sleep to make sure we don't ping slurm too frequently
	os.system("sleep 1")

	#read the slurm queue
	slurm_queue = open("slurm_queue.txt", "r")
	for slurm_line in slurm_queue.readlines():
		#extract the first entry in the line (which would be the id), and determine if it is in the list of active slurm jobs
		#if it is, set jobs running to true, because this id is in the list of active jobs
		cur_id = slurm_line.split()[0]

		if cur_id in active_slurm_jobs:
			jobs_running = True