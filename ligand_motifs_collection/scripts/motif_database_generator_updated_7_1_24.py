#the goal of this program is to take a directory of ligand-containing protein pdb files, and extract motifs
#the program needs to generate param files of the ligands, and will then run the ligand_motifs.cc Rosetta script to build a motifs library
#motif database will be created in the same directory as the pdb database
#this program and its supporting auto_params files are to be in the directory of the pdb  database

#run the program as: python motif_database_generator.py 

#import packages:
import sys
import os
import argparse

# Create the parser
parser = argparse.ArgumentParser(
	description="This script utilizes a Slurm job system to automate the collection of motifs from a directory of protein-ligand pdb files from RCSB. This script requires a path to "
)

# Add the rosetta executable path argument
parser.add_argument(
	'-r', '--rpath',
	type=str,
	required=True,
	help='Full path to the Rosetta executable for ligand motif collection. Within Rosetta, the path is /rosetta/source/bin/identify_ligand_motifs.linuxgccrelease.'
)

# Add the auto_params path argument
parser.add_argument(
	'-a', '--apath',
	type=str,
	required=True,
	help='Full path to the auto_params.sh executable.'
)

# Add the SeparatePDBsByChain path argument
parser.add_argument(
	'-s', '--spath',
	type=str,
	required=True,
	help='Full path to the SeparatePDBsByChain.pl executable.'
)

# Add the working location path argument
parser.add_argument(
	'-w', '--wpath',
	type=str,
	required=False,
	help='Full path to the working location where the pdb files are. This is optional, and will use the current location if not specified.'
)

# Add the working location path argument
parser.add_argument(
	'-m', '--mpath',
	type=str,
	required=False,
	help='Full path to the location where the Rosetta molfile_to_params.py script is.'
)

#parse arguments
args = parser.parse_args()

#set arguments to variables and add a backslash to the end if one was not included
rosetta_path = args.rpath
if rosetta_path.endswith("/") == False:
	rosetta_path = rosetta_path + "/"
auto_path = args.apath
if auto_path.endswith("/") == False:
	auto_path = auto_path + "/"
separate_path = args.spath
if separate_path.endswith("/") == False:
	separate_path = separate_path + "/"
mol_to_param_path = args.mpath
if mol_to_param_path.endswith("/") == False:
	mol_to_param_path = mol_to_param_path + "/"

#initial setup of working path
working_path = ""

#set up location variable to be either working path or current location,
location = os.getcwd()

#set up working path if -w or --wpath was used
if '-w' in vars(args) or '--wpath' in vars(args):
	working_path = args.spath

	#remove backslash from end of working path if it is there
	if working_path.endswith("/") == True:
		working_path = working_path[:-1]

	location = working_path

	#move to the location
	os.chdir(location)


#identify user name via whoami for use with slurm queue browsing
#user_id = os.getlogin()
user_id = os.path.expanduser('~').split("/")[len(os.path.expanduser('~').split("/")) - 1]

#look for all pdb files in directory, compose a list of the pdbs and run auto
pdb_list = []

for r, d, f in os.walk(location):
	for file in f:
		if file.endswith(".pdb"):
			#append the pdb name, which is whatever comes before the period (in case there are pdbs that are more than just their 4 character code)
			pdb_list.append(file.split(".")[0])

#run auto_params.sh for each pdb file in the pdb list
#put the generated files in their own folder by pdb
#iterate the loop to perform this operation for groups of 100 files

"""
#old method
for file in pdb_list:
	#first 4 characters are the pdb to make the folder
	os.system("mkdir " + file)

	#run auto_params
	os.system("./auto_params.sh " + file + ".pdb")
	
	#move all files for the pdb to the folder
	os.system("mv " + file + "* " + file)
"""

#list to hold up to 100 pdbs from the list
mini_list = []

for index in range(len(pdb_list)):
	#add the pdb from pdb_list[index] to the mini_list of up to 100 pdbs
	#only add if the pdb does not contain "Ligatoms" or is a "pdb_pdb" (i.e. 1a1f_1a1f.pdb) name
	mini_list.append(pdb_list[index])

	#when there are 100 entries in mini_list, make a job file and clear mini_list to repeat
	if (len(mini_list) % 100) == 0:
		#write the job file that calls auto_params
		auto_params_job_file_name = location + "/" + str(len(mini_list)) + ".job"
		auto_params_job_writer = open(auto_params_job_file_name, "w")
		#write header lines for job
		auto_params_job_writer.write("#!/bin/bash\n")
		auto_params_job_writer.write("#SBATCH -p short # Partition to submit to\\n\")\n")
		auto_params_job_writer.write("#SBATCH -n 1 # Number of cores requested\\n\")\n")
		auto_params_job_writer.write("#SBATCH -N 1 # Ensure that all cores are on one machine\\n\")\n")
		auto_params_job_writer.write("#SBATCH -t 720 # Runtime in minutes\\n\")\n")
		auto_params_job_writer.write("#SBATCH --mem=12800 # Memory per cpu in MB (see also --mem-per-cpu)\\n\")\n")
		auto_params_job_writer.write("#SBATCH -o " + location + "/" + str(index) + "_hostname_%A_%a.out # Standard out goes to this file\\n\")\n")
		auto_params_job_writer.write("#SBATCH -e " + location + "/" + str(index) + "_hostname_%A_%a.err # Standard err goes to this filehostname\\n\")\n")
		auto_params_job_writer.write("cd " + location + "/" + "\n")

		#write auto_params call for each pdb in mini_list
		for item in mini_list:
			#first 4 characters are the pdb to make the folder
			auto_params_job_writer.write("mkdir " + item + "\n")
			#run auto_params
			auto_params_job_writer.write(auto_path + "/auto_params.sh " + item + ".pdb " + separate_path + " " + mol_to_param_path + "\n")
			#move all files for the pdb to the folder
			auto_params_job_writer.write("mv " + item + "* " + item + "\n")

		auto_params_job_writer.close()

		#clear contents of mini_list
		mini_list = []

		#send off the generated job file
		os.system("sbatch " + auto_params_job_file_name)

#send off one last call of auto_params for remaining contents of mini_list
auto_params_job_file_name = location + "/" + str(len(mini_list)) + ".job"
auto_params_job_writer = open(auto_params_job_file_name, "w")
#write header lines for job
auto_params_job_writer.write("#!/bin/bash\n")
auto_params_job_writer.write("#SBATCH -p short # Partition to submit to\\n\")\n")
auto_params_job_writer.write("#SBATCH -n 1 # Number of cores requested\\n\")\n")
auto_params_job_writer.write("#SBATCH -N 1 # Ensure that all cores are on one machine\\n\")\n")
auto_params_job_writer.write("#SBATCH -t 720 # Runtime in minutes\\n\")\n")
auto_params_job_writer.write("#SBATCH --mem=1280 # Memory per cpu in MB (see also --mem-per-cpu)\\n\")\n")
auto_params_job_writer.write("#SBATCH -o " + location + "/" + "end" + "_hostname_%A_%a.out # Standard out goes to this file\\n\")\n")
auto_params_job_writer.write("#SBATCH -e " + location + "/" + "end" + "_hostname_%A_%a.err # Standard err goes to this filehostname\\n\")\n")
auto_params_job_writer.write("cd " + location + "/" + "\n")

#write auto_params call for each pdb in mini_list
for item in mini_list:
	#first 4 characters are the pdb to make the folder
	auto_params_job_writer.write("mkdir " + item + "\n")
	#run auto_params
	auto_params_job_writer.write(auto_path + "/auto_params.sh " + item + ".pdb " + separate_path + " " + mol_to_param_path + "\n")
	#move all files for the pdb to the folder
	auto_params_job_writer.write("mv " + item + "* " + item + "\n")

auto_params_job_writer.close()

#clear contents of mini_list
mini_list = []

#send off the generated job file
os.system("sbatch " + auto_params_job_file_name)

#stall the program until all auto_params calls are done
#accomplished by stalling until all pdb files are absent from the working directory

#variable to indicate whether pdbs are present in the working directory (diagnostic of auto_params not being done)
pdbs_present = True

while pdbs_present == True:
	#set value of pdbs_present to False, will turn true if pdbs present
	pdbs_present = False

	#get the contents of the directories
	for r, d, f in os.walk(location):
		

		#look only at the root directory that is the working directory
		if r == location:
			#look through all files to see if there is a .pdb file
			for file in f:
				if file.endswith(".pdb"):
					#print(file)
					pdbs_present = True

					#remove the rogue "_.pdb" file that hangs up the program
					if file == "_.pdb":
						os.system("mkdir _")
						os.system("mv _.pdb _")

#Build a slurm job and flag file for each param/pdb file generated
for pdb in pdb_list:

	#make list of param files
	param_list = []

	#walk the pdb directory to get all params files
	for r, d, f in os.walk(location + "/" + pdb):
		for file in f:
			if file.endswith(".params"):
				#append the params name, which is whatever comes before the period (in case there are pdbs that are more than just their 4 character code)
				param_list.append(file.split(".")[0])

	#make a pdb file for each ligand (params file), where there is only 1 ligand in the pdb file
	#-------split pdb
	#read in the base pdb file
	pdb_file = open(location + "/" + pdb + "/" + pdb + ".pdb", "r")

	#blank list of ligands
	ligands = []

	#read through the file and get the list of ligands
	for line in pdb_file.readlines():
		if line.startswith("HETATM"):

			#get the ligand name
			lig = line[17:20].strip(" ")
			
			#ignore water
			if  lig != "HOH":
				#check the list of ligand, and if lig is not in the list, add it to the list

				if lig not in ligands:
					ligands.append(lig)

	#print(ligands)
	pdb_file.close()

	#empty list to hold single ligand pdb files
	single_ligand_files  = []

	#make a pdb,  ignoring all except the focused ligand
	for ligand in ligands:
		#open a file to write the new pdb to
		write_file = open(location + "/" + pdb + "/" + pdb + "_" + ligand +  ".pdb", "w")
		single_ligand_files.append(location + "/" + pdb + "/" + pdb + "_" + ligand +  ".pdb")
		pdb_file = open(location + "/" + pdb + "/" + pdb + ".pdb", "r")
		#run  through each line of the read pdb, and only write lines that do not contain the other ligands
		for line in pdb_file.readlines():

			write_line = True

			#check if the line has any ligand in it other than the desired ligand
			for lig in ligands:
				if lig != ligand and lig in line:
					write_line = False

			if write_line == True:
				#print(line)
				write_file.write(line)

		write_file.close()
		pdb_file.close()
	#---end split

	#make a flag file for each split pdb file file
	for file in single_ligand_files:
		#"file_flags is the flag file name"

		file_name = file.split("/")[len(file.split("/")) - 1].split(".")[0]

		flag_file_name = location + "/" + pdb + "/" + file_name + "_flags"
		flag_writer = open(flag_file_name, "w")
		#write contents of file
		#flag_writer.write("-extra_res_fa " + location + "/" + pdb + "/" + file  + ".params\n")
		#flag_writer.write("-extra_res_fa ")

		for param_file in param_list:
			flag_writer.write(param_file  + ".params ")
		
		flag_writer.write("\n")
		flag_writer.write("-ignore_unrecognized_res\n")
		flag_writer.write("-hb_score_cutoff -0.3\n")
		flag_writer.write("-pack_score_cutoff -0.5\n")
		#flag_writer.write("-s " + location + "/" + pdb + "/" + pdb + ".pdb\n")
		split_file = file.split("/")[len(file.split("/")) - 1]
		flag_writer.write("-s " + split_file + "\n")
		flag_writer.write("-output_motifs_as_pdb false \n")
		flag_writer.close()

	#make and submit a slurm job for each collective of flag files
	job_file_name = location + "/" + pdb + "/" + pdb + ".job"
	job_writer = open(job_file_name, "w")

	#write header lines for job
	job_writer.write("#!/bin/bash\n")
	job_writer.write("#SBATCH -p short # Partition to submit to\\n\")\n")
	job_writer.write("#SBATCH -n 1 # Number of cores requested\\n\")\n")
	job_writer.write("#SBATCH -N 1 # Ensure that all cores are on one machine\\n\")\n")
	job_writer.write("#SBATCH -t 720 # Runtime in minutes\\n\")\n")
	job_writer.write("#SBATCH --mem=1280 # Memory per cpu in MB (see also --mem-per-cpu)\\n\")\n")
	job_writer.write("#SBATCH -o " + location + "/" + pdb + "/" + pdb + "_hostname_%A_%a.out # Standard out goes to this file\\n\")\n")
	job_writer.write("#SBATCH -e " + location + "/" + pdb + "/" + pdb + "_hostname_%A_%a.err # Standard err goes to this filehostname\\n\")\n")
	job_writer.write("cd " + location + "/" + pdb + "\n")

	#add lines to the file to run each flag file, and handle the aftermath of each run before running the next
	for file in single_ligand_files:
		file_name = file.split("/")[len(file.split("/")) - 1].split(".")[0]

		flag_file_name = location + "/" + pdb + "/" + file_name + "_flags"
		job_writer.write("mkdir " + location + "/" + pdb + "/Ligand_motif_dir\n")
		#run ligand_motifs
		job_writer.write(rosetta_path + "identify_ligand_motifs.linuxgccrelease @" + flag_file_name + "\n")
		#make the "Ligand_motif_dir" the name of the param file
		job_writer.write("mv Ligand_motif_dir " + file_name + "\n")
		#move AllMattMotifsFile.motifs to the renamed folder
		job_writer.write("mv AllMattMotifs.motifs " + file_name + "\n")
	

	job_writer.close()


	#check queue length to see how many jobs are in queue, do not submit job until queue length is 400 or less
    #make file to hold length of squeue
    
	os.system("squeue -A " + user_id + " | wc -l  > squeue_file.txt")
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


	        os.system("squeue -A " + user_id + " | wc -l  > squeue_file.txt")
	        squeue_length_file = open("squeue_file.txt", 'r')
	        squeue_length = ""
	        for line in squeue_length_file.readlines():
	                #rebuild squeue_length to ensure it is only numbers
	                for char in line:
	                        if char.isnumeric():
	                                squeue_length = squeue_length + char
	        squeue_length = int(squeue_length)
	        squeue_length_file.close()

	#submit the job
	os.system("sbatch " + job_file_name)
