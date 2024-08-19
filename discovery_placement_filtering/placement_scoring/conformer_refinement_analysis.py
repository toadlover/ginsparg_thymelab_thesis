#this script looks at a a csv file of placement score data (i.e. raw_scores.csv) to observe a list of placement pdbs and compare them to sets of placements from expanded conformers to determine if hte expanded conformers improve the initial placement
#this script requires a raw_scores.csv (or similar layout file) and a directory to where raw_scores.csv files are for the conformer placements of the initial ligand

#i.e. for the ligand PV-006690961802 in the input raw_scores.csv file, this script will look in the provided path for another raw_scores.csv file for PV-006690961802 conformers that were placed, and compare the placements
#comparisons include the following (with a difference value calculated for each):
#rmsd from initial placements (ensure that we are placing similarly)
#Rosetta ddg
#total motif-like interactions found
#total real motifs found
#total Rosetta hbonds found
#total Rosetta hbond energy
#closest autodock recovery
#closest autodock recovery ddg
#torsion strain energy

#for each comparison, determine which are improvements and get a total of the number of metrics improved (represented as ratio of improved over total) for each placement that can be used as a quick reference
#list metrics for all placements and highlight all that improve (more help that hurt)

import os,sys

#get the initial csv
placements_csv = open(sys.argv[1], "r")

#get the location where palcements are
conf_placements_locations = sys.argv[2]

#add a slash to the location if one does not exist
if conf_placements_locations.endswith("/") == False:
	conf_placements_locations = conf_placements_locations + "/"

#move to the conf_placement_locations
os.chdir(conf_placements_locations)

#make a file to hold the conf delta data and data on the count of placements for each conformer
delta_data_file = open("delta_scores.csv","w")
#delta file with the most improved placement
delta_data_file_most_improved = open("delta_scores_most_improved.csv","w")
#placements_per_initial_conf_file = open("placements_per_initial_conformer.csv", "w")

#list to hold the score term names
score_term_names = []
score_term_names_delta_write_file = []

#read placements csv
for line in placements_csv.readlines():
	#handle the header and create score term lines
	if line.startswith("file,"):
		for item in line.split(","):
			score_term_names.append(item.strip())


		#write header lines for the delta data file
		#delta_data_file.write("initial_file,compare_file,initial_ddg,compare_ddg,delta_ddg,initial_total_motifs,compare_total_motifs,delta_total_motifs,initial_significant_motifs,compare_significant_motifs,delta_significant_motifs,initial_real_motif_ratio,hbond_motif_count,hbond_motif_energy_sum,closest_autodock_recovery_rmsd,closest_autodock_recovery_ddg,strain_energy,total")
		for item in score_term_names:
			#add a value for the initial placement, the ocmpare from the conformer, and if not one of values that do not involve delta, the delta comparison between values
			#nondelta values are the filenames, rmsd, and the count for the number of improved values
			score_term_names_delta_write_file.append(item + "_initial")
			score_term_names_delta_write_file.append(item + "_compare")
			if item != "file":
				score_term_names_delta_write_file.append(item + "_delta")

		#add values for rmsd and improvement ratio
		score_term_names_delta_write_file.append("rmsd")
		score_term_names_delta_write_file.append("improvement_ratio")

		#write the terms to the delta files
		for item in score_term_names_delta_write_file:
			delta_data_file.write(item + ",")
			delta_data_file_most_improved.write(item + ",")

		delta_data_file.write("\n")
		delta_data_file_most_improved.write("\n")

		continue

	#handle a new conformer from the main list
	#create a dictionary for this ligand with its corresponding values and fill with values
	initial_conf_dict = {}

	#counter for items to apply the name to the score term in the dictionary
	term_counter = 0

	for item in line.split(","):
		score_term_names.append(item.strip())
		#add the term to the dictionary
		initial_conf_dict[score_term_names[term_counter]] = item.strip()

		term_counter = term_counter + 1

		#extract the ligand name
	#ligand name is the 3rd from last when splitting by underscores
	lig_name = initial_conf_dict["file"].split("_")[len(initial_conf_dict["file"].split("_")) - 3]
	print(lig_name)

	#get the coordinates of the atoms in the initial conformer for an rmsd comparison
	#ignore hydrogens
	#Store in dictionary with atoms as keys
	initial_conf_atom_coords = {}

	#open the pdb, which is under the file keyname
	initial_conf_pdb = open(initial_conf_dict["file"],"r")

	for line2 in initial_conf_pdb.readlines():
		#HETATM lines have the ligand data
		if line2.startswith("HETATM "):
			#3rd position (index value 2 by splitting)
			unique_atom = line2.split()[2]

			#continue if Hydrogen
			#easy solution is if the first character is H, since n oother relevant element starts with H
			if unique_atom[0] == "H":
				continue

			#otherwise add the atom to the dictionary with xyz values as list
			initial_conf_atom_coords[unique_atom] = [float(line2.split()[6]),float(line2.split()[7]),float(line2.split()[8])]

	initial_conf_pdb.close()

	"""
	#with data read in, now find the raw_scores csv that coorresponds to the conformer and read it to perform the comparison analysis
	for r,d,f in os.walk(conf_placements_locations):
		#iterate over files
		for file in f:
			#if the ligand name is in the root and the file is raw_scores.csv
			if lig_name in r and file == "raw_scores.csv":
				#read the file so we can compare
	"""

	#set up a variable to hold the best placement
	best_placement = ""

	#compare_file = open(r + "/raw_scores.csv", "r")
	compare_file = open(conf_placements_locations + lig_name + "/placements/raw_scores.csv", "r")

	for compare_line in compare_file.readlines():
		
		#skip the header line since the header should be the same as the initial file
		if compare_line.startswith("file"):					
			continue

		#set up dictionary for the score terms from the compare file
		compare_conf_dict = {}
		term_counter = 0
		for item in line.split(","):
			score_term.append(item.strip())
			#add the term to the dictionary
			compare_conf_dict[score_term_names[term_counter]] = item.strip()

			term_counter = term_counter + 1

		#get the delta between like terms (besides file) for the initial and compare, and determine if the compare is an improvement
		term_counter = 0
		improved_term_counter = 0

		for i in range(len(score_term_names_delta_write_file)):
			#write the comparisons
			#handle "file" differently since it is a string
			if "file" in score_term_names_delta_write_file[i]:
				#write the relevant file
				if "initial" in score_term_names_delta_write_file[i]:
					delta_data_file.write(initial_conf_dict["file"] + ",")
				if "compare" in score_term_names_delta_write_file[i]:
					delta_data_file.write(compare_conf_dict["file"] + ",")

				continue

			#handle rmsd
			if score_term_names_delta_write_file[i] == "rmsd":
				#hard code to 0 for now for testing, will fix later
				delta_data_file.write("0,")
				continue

			#handle the improvement ratio
			if score_term_names_delta_write_file[i] == "improvement_ratio":
				#hard code to 0 for now, and we will implement later
				delta_data_file.write("0,")
				continue

			#handle the other comparison terms

			#initial/compare/delta off the term for accessing in the corresponding dictionary
			stripped_term = score_term_names_delta_write_file[i].split("_initial")[0].split("_compare")[0].split("_delta")[0]

			if "initial" in score_term_names_delta_write_file[i]:
				delta_data_file.write(str(initial_conf_dict[stripped_term]) + ",")
			if "compare" in score_term_names_delta_write_file[i]:
				delta_data_file.write(str(compare_conf_dict[stripped_term]) + ",")
			#if delta
			if "delta" in score_term_names_delta_write_file[i]:

				#delta is initial minus compare
				delta_value = float(initial_conf_dict[stripped_term]) - float(compare_conf_dict[stripped_term])

				delta_data_file.write(str(delta_value) + ",")

		#write a newline
		delta_data_file.write("\n")









