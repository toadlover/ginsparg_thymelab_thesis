#script to look at a directory of pdbs of placed ligands using Rosetta's liganddiscoverysearch
#ligands are scored based on values derived from rosetta (ddg, motif counts and real ratio), ligand strain energy with tldr, and attepmter recovery with AutoDock vina
#system scores are summative of all scoring values used to calculate a score, with a larger and more positive value ideally being more favorable

#package imports
import os,sys
import argparse

# Create the parser
parser = argparse.ArgumentParser(
	description="This script looks at a directory of placed ligand pdb files from Rosetta's liganddiscovery search. It scores ligands based on values from Rosetta (ddg, motif counts and real ratio), ligand strain energy using tldr STRAIN, and attempted recovery of the Rosetta prediction with AutoDock Vina."
)

# Add the working location path argument
#optional; if not used, will use current location to look for pdb files
parser.add_argument(
	'-w', '--wpath',
	type=str,
	required=False,
	help='(Optional) Full path to the working location where the pdb files are. Will use the current location if not specified.'
)

# Add the autodock automation script path argument
#run_autodock_on_placed_ligands.py
#optional, will not attempt autodock if flag is not used; also requires the use of the -v flag that indicates where the build of autodock vina is
parser.add_argument(
	'-a', '--apath',
	type=str,
	required=False,
	help='(Optional) Full path to the location where the run_autodock_on_placed_ligands.py script is located for attempting to recover Rosetta placements with AutoDock Vina. Recovery will not be attempted if this flag is not used or if the path to the AutoDock Vina Executable (-v, --vpath) is not used. This is designed and tested on Linux using AutoDock Vina version 1.1.2, and may not work on other platforms.'
)

# Add the AutoDock Vina executable path argument
#looking for executable called "vina"
#optional, will not attempt autodock if flag is not used; also requires the use of the -a flag that indicates where the autodock automation script is
parser.add_argument(
	'-v', '--vpath',
	type=str,
	required=False,
	help='(Optional) Full path to the location where the AutoDock Vina "vina" executable is located for attempting to recover Rosetta placements with AutoDock Vina. Recovery will not be attempted if this flag is not used or if the path to the AutoDock Vina automation script run_autodock_on_placed_ligands.py (-a, --apath) is not used. This is designed and tested on Linux using AutoDock Vina version 1.1.2, and may not work on other platforms or versions.'
)

# Add the tldr STRAIN executable path argument
#looking for a script called run_torsion_check_on_placed_ligands.py
#optional, will not attempt strain if  if flag is not used; also requires the use of the -t flag that indicates where the autodock automation script is
parser.add_argument(
	'-r', '--rpath',
	type=str,
	required=False,
	help='(Optional) Full path to the location where the run_torsion_check_on_placed_ligands.py script is located for evaluating conformer strain energies. Will not be attempted if this flag is not used or if the path to the tldr strain automation script run_autodock_on_placed_ligands.py (-t, --tpath) is not used. This is designed and tested on Linux, and may not work on other platforms. The script requires a python or conda environment to have openbabel (obabel) to be installed.'
)

# Add the tldr STRAIN executable path argument
#looking for a script called Torsion_Strain.py
#optional, will not attempt strain if  if flag is not used; also requires the use of the -r flag that indicates where the autodock automation script is
parser.add_argument(
	'-t', '--tpath',
	type=str,
	required=False,
	help='(Optional) Full path to the location where the tldr strain Torsion_Strain.py script is located for evaluating conformer strain energies. Will not be attempted if this flag is not used or if the path to the tldr strain automation script run_autodock_on_placed_ligands.py (-r, --rpath) is not used. This is designed and tested on Linux, and may not work on other platforms. The script requires a python or conda environment to have rdkit to be installed.'
)

# Add the score_weights file path argument
#looking for a file called score_weights.csv
#optional, collects user-set score weights for values used on scoring placed ligands; any value that does not have a weight set is defaulted to 1
parser.add_argument(
	'-s', '--weights_path',
	type=str,
	required=False,
	help='(Optional) Full path to the location where the score weights file score_weights.csv is located. The format within the file is that each line has a score term at index 0 and the score weight at index 1.'
)

#parse arguments
args = parser.parse_args()

#debug print of args
print(args)

#initial setup of working path
working_path = ""

#set up location variable to be either working path or current location,
location = os.getcwd()

#set up working path if -w or --wpath was used
#if 'w' in vars(args) or 'wpath' in vars(args):
if args.wpath != None:
	working_path = args.wpath

	#remove backslash from end of working path if it is there
	if working_path.endswith("/") == True:
		working_path = working_path[:-1]

	location = working_path

	#move to the location
	os.chdir(location)

#determine whether to use autodock and set up autodock variables
autodock_exec_path = ""
autodock_automate_path = ""
#default true, set to false if either relevant flag is not used
attempt_autodock = True
#autodock automate path
#if 'v' in vars(args) or 'vpath' in vars(args):
if args.vpath != None:
	autodock_exec_path = args.vpath
else:
	attempt_autodock = False

#autodock automate path
#if 'a' in vars(args) or 'apath' in vars(args):
if args.apath != None:
	autodock_automate_path = args.apath
	#print(autodock_automate_path)
	#add backslash to end of path if there is not one
	if autodock_automate_path.endswith("/") == False:
		autodock_automate_path = autodock_automate_path + "/"
	#print(autodock_automate_path)
else:
	attempt_autodock = False

#print('-v' in vars(args), '--vpath' in vars(args), '-v' in vars(args) or '--vpath' in vars(args), '-a' in vars(args), '--apath' in vars(args), '-a' in vars(args) or '--apath' in vars(args))
#print('v' in vars(args), 'vpath' in vars(args), 'v' in vars(args) or 'vpath' in vars(args), 'a' in vars(args), 'apath' in vars(args), 'a' in vars(args) or 'apath' in vars(args))

#determine whether to use tldr strain and set up variables
strain_exec_path = ""
strain_automate_path = ""
#default true, set to false if either relevant flag is not used
attempt_strain = True
#strain torsion path
#if 't' in vars(args) or 'tpath' in vars(args):
if args.tpath != None:
	strain_exec_path = args.tpath

	#add backslash to end of path if there is not one
	if strain_exec_path.endswith("/") == False:
		strain_exec_path = strain_exec_path + "/"
else:
	attempt_strain = False

#strain automate path
#if 'r' in vars(args) or 'rpath' in vars(args):
if args.rpath != None:
	strain_automate_path = args.apath

	#add backslash to end of path if there is not one
	if strain_automate_path.endswith("/") == False:
		strain_automate_path = strain_automate_path + "/"
else:
	attempt_strain = False

#prepare score weights, set up initial score weight dictionary
#if adjusting weights in score_weights.csv, any changed weight must match one of these entries exactly (including case)
#if a value's weight is set to zero and would otherwise be used, it will not be considered in scoring
score_weights = {"ddg":1,"total_motifs":1,"significant_motifs":1,"real_motif_ratio":1,"closest_autodock_recovery_rmsd":1,"closest_autodock_recovery_ddg":1,"strain_energy":1}

#if weights_path flag is set, work on adjusting weights
#if 's' in vars(args) or 'weights_path' in vars(args):
if args.weights_path != None:
	weights_path = args.weights_path

	#add backslash to end of path if there is not one
	if weights_path.endswith("/") == False:
		weights_path = weights_path + "/"

	#open the weights file and read the weights data
	weights_file = open(weights_path + "score_weights.csv", "r")

	#read each line and get the score term and weight
	#when splitting line by "," the score term should be index 0 and score be index 1
	#not all score terms need to be in the csv file, and order of terms does not matter
	#there is no header line
	#i.e.:
	"""
	ddg,3
	real_motif_ratio,10
	closest_autodock_recovery_ddg,0
	strain_energy,0.5
	"""
	for line in weights_file.readlines():
		term = line.split(",")[0]
		weight = float(line.split(",")[1].strip())
		
		#adjust the term
		if term in score_weights.keys():
			score_weights[term] = weight
		else:
			print("WARNING: included term '" + term + "' is not valid!")
			print("Valid terms are ddg, total_motifs, significant_motifs, real_motif_ratio, closest_autodock_recovery_rmsd, closest_autodock_recovery_ddg, and strain_energy.")

#set up a file to hold the scoring for the placed ligands in csv format both with weighted and non-weighted values
#write the file in the location
#raw score file, generally for debugging and score term weights are not applied
raw_score_file = open(location + "/raw_scores.csv", "w")
#write header line
raw_score_file.write("file,ddg,total_motifs,significant_motifs,real_motif_ratio,closest_autodock_recovery_rmsd,closest_autodock_recovery_ddg,strain_energy,total\n")
#set up weighted score file; this file should be used for selection
weighted_score_file = open(location + "/weighted_scores.csv", "w")
#write header line
weighted_score_file.write("file,ddg,total_motifs,significant_motifs,real_motif_ratio,closest_autodock_recovery_rmsd,closest_autodock_recovery_ddg,strain_energy,total\n")

#look at all pdb files in the given directory and get started on scoring
for r,d,f in os.walk(location):
	for file in f:
		#only look in the current location (do not look further down)
		if r != location:
			continue

		#only look at pdb files, and make sure to ignore the minipose file
		if file.endswith(".pdb") and "minipose" not in file:
			#set initial values for score terms when scored and unscored
			#set as 2 entry list initially seeded as ""
			#entry 0 is the raw value and entry 1 is with the weight applied
			ddg = [0,0]
			total_motifs = [0,0]
			significant_motifs = [0,0]
			real_motif_ratio = [0,0]
			#raw starting value set different for selecting closest rmsd value
			closest_autodock_recovery_rmsd = [100,0]
			closest_autodock_recovery_ddg = [0,0]
			strain_energy = [0,0]
			total = [0,0]

			#make a directory that contains the starting pdb file and will hold any intermediate data
			dir_name = file.split(".")[0]
			os.system("mkdir " + dir_name)
			os.system("cp " + file + " " + dir_name)

			#extract comment table data from the file to potentially assign values for ddg, total_motifs, significant_motifs, and real_motif_ratio
			#read the file
			placement_file = open(r + "/" + file, "r")

			#get the data
			for line in placement_file.readlines():
				#ddg
				if line.startswith("Scoring: Post-HighResDock system ddG:"):
					ddg = [float(line.split()[len(line.split()) - 1].strip()),float(line.split()[len(line.split()) - 1].strip())*score_weights["ddg"]]
				#total_motifs
				if line.startswith("Placement motifs: Total motifs made:"):
					total_motifs = [float(line.split()[len(line.split()) - 1].strip()),float(line.split()[len(line.split()) - 1].strip())*score_weights["total_motifs"]]
				#significant_motifs
				if line.startswith("Placement motifs: Motifs made against significant residues count:"):
					significant_motifs = [float(line.split()[len(line.split()) - 1].strip()),float(line.split()[len(line.split()) - 1].strip())*score_weights["significant_motifs"]]
				#real_motif_ratio
				if line.startswith("Placement motifs: Real motif ratio:"):
					real_motif_ratio = [float(line.split()[len(line.split()) - 1].strip()),float(line.split()[len(line.split()) - 1].strip())*score_weights["real_motif_ratio"]]

			#attempt autodock if able
			if attempt_autodock:

				print("Attempting autodock")

				#call run_autodock_on_placed_ligands.py on the directory for this file
				#path to automation script, path to working file, path to autodock executable
				os.system("python " + autodock_automate_path + "run_autodock_on_placed_ligands.py " + r + "/" + dir_name + " " + autodock_exec_path)

				#read the corresponding autodock_data.csv file and determine the closest placement if there is any (and get corresponding energy score)
				#if there are no placement attempts, set values for closest_autodock_recovery_rmsd to 100 (arbitrary high value that should tank the placement because getting no placement attempts is bad) and closest_autodock_recovery_ddg to 0
				autodock_file = open(dir_name + "/autodock_data.csv", "r")

				for line in autodock_file.readlines():
					#ignore the header line
					if "system,model,rmsd,close_to_rosetta,vina_energy" in line:
						continue
					#handle if there were no placements; if split(,)[1] is 0
					if line.split(",")[1] == "0":
						closest_autodock_recovery_ddg = [0,0*score_weights["closest_autodock_recovery_ddg"]]
						closest_autodock_recovery_rmsd = [100,100*score_weights["closest_autodock_recovery_rmsd"]]
					#otherwise, handle the data extract the closest rmsd with its corresponding energy
					if float(line.split(",")[2]) < closest_autodock_recovery_rmsd[0]:
						#set the new closest rmsd with its ddg
						closest_autodock_recovery_rmsd = [float(line.split(",")[2]), float(line.split(",")[2]) * score_weights["closest_autodock_recovery_rmsd"]]
						closest_autodock_recovery_ddg = [float(line.split(",")[4].strip()), float(line.split(",")[4].strip()) * score_weights["closest_autodock_recovery_ddg"]]

			#attempt strain torsion if able
			if attempt_strain:
				
				print("Attempting tldr STRAIN")

				#call the run_torsion_check_on_placed_ligands.py script
				print("python " + strain_automate_path + "run_torsion_check_on_placed_ligands.py " + r + "/" + dir_name + " " + strain_exec_path)
				os.system("python " + strain_automate_path + "run_torsion_check_on_placed_ligands.py " + r + "/" + dir_name + " " + strain_exec_path)

				#read the torsion data table csv
				#(pdb name)_torsion_data.csv
				strain_file = open(dir_name + "/" + dir_name + "_lig_Torsion_Strain.csv", "r")

				#read the strain file and get the total strain energy if it was calculated
				for line in strain_file.readlines():
					#there is a chance the the energy is inexplicably not calculated
					#behavior for this case
					#split(",")[1] of line is "NA" in this case
					if line.split(",")[1].strip() == "NA":
						#set strain energy to 100
						strain_energy = [100,100*score_weights["strain_energy"]]
					else:
						#get average of high and low energy predictions to use
						high = float(line.split(",")[1])
						low = float(line.split(",")[2])
						average = (high + low) / 2

						strain_energy = [average,average*score_weights["strain_energy"]]

			total[0] = ddg[0] + total_motifs[0] + significant_motifs[0] + real_motif_ratio[0] + closest_autodock_recovery_rmsd[0] + closest_autodock_recovery_ddg[0] + strain_energy[0]
			total[1] = ddg[1] + total_motifs[1] + significant_motifs[1] + real_motif_ratio[1] + closest_autodock_recovery_rmsd[1] + closest_autodock_recovery_ddg[1] + strain_energy[1]

			#write scores to the score files
			#file,ddg,total_motifs,significant_motifs,real_motif_ratio,closest_autodock_recovery_rmsd,closest_autodock_recovery_ddg,strain_energy,total
			raw_score_file.write(file + "," + str(ddg[0]) + "," + str(total_motifs[0]) + "," + str(significant_motifs[0]) + "," + str(real_motif_ratio[0]) + "," + str(closest_autodock_recovery_rmsd[0]) + "," + str(closest_autodock_recovery_ddg[0]) + "," + str(strain_energy[0]) + "," + str(total[0]) + "\n")
			weighted_score_file.write(file + "," + str(ddg[1]) + "," + str(total_motifs[1]) + "," + str(significant_motifs[1]) + "," + str(real_motif_ratio[1]) + "," + str(closest_autodock_recovery_rmsd[1]) + "," + str(closest_autodock_recovery_ddg[1]) + "," + str(strain_energy[1]) + "," + str(total[1]) + "\n")
