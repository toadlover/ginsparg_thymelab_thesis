#this is a modified version of conformer_refinement_analysis that is meant to operate  on a single system as opposed to a directory of systems. This can be used for parallelization (such as a slurm array job)
#this allows you to work more broadly and instead of a folder of refined placements, you just run the script in 2 placement csv files to compare against the two

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

#commenting old code for in input location
"""
#get the location where palcements are
conf_placements_locations = sys.argv[2]

#add a slash to the location if one does not exist
if conf_placements_locations.endswith("/") == False:
	conf_placements_locations = conf_placements_locations + "/"

#move to the conf_placement_locations
os.chdir(conf_placements_locations)
"""
#instead we take in a second placements csv file (such as a raw_scores.csv)
#can include a path
refined_csv = sys.argv[2]

#make a file to hold the conf delta data and data on the count of placements for each conformer
delta_data_file = open("delta_scores.csv","w")
#delta file with the most improved placement
delta_data_file_most_improved = open("delta_scores_most_improved.csv","w")
#placements_per_initial_conf_file = open("placements_per_initial_conformer.csv", "w")

#list to hold the score term names
score_term_names = []
score_term_names_delta_write_file = []

#placement counter for tracking
placement_counter = 0

#read placements csv
for line in placements_csv.readlines():
	
	placement_counter = placement_counter + 1

	#handle the header and create score term lines
	if line.startswith("file,"):
		for item in line.split(","):
			score_term_names.append(item.strip())


		#if we have a ligand column, we can just directly add it to the delta terms and not do initial or compare, since it should come out the same
		if item == "ligand":
			score_term_names_delta_write_file.append(item)

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
	print(lig_name, placement_counter)

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

	#set up a list to hold the best placement(s) and all placements
	#best_placements = []

	all_placements = []

	#hold the best improvement ratio value to help select which placements are the best
	best_improvement_ratio = 0

	#sanity check to make sure that the raw scores file exists
	"""
	if os.path.exists(conf_placements_locations + lig_name + "/placements/raw_scores.csv") == False:
		print("raw score file does not exist...")
		continue

	#compare_file = open(r + "/raw_scores.csv", "r")
	compare_file = open(conf_placements_locations + lig_name + "/placements/raw_scores.csv", "r")
	"""

	compare_file = open(refined_csv , "r") 

	for compare_line in compare_file.readlines():
		
		#skip the header line since the header should be the same as the initial file
		if compare_line.startswith("file"):					
			continue

		#ensure that we are only comparing the placements of the same ligand
		#break,
		if lig_name not in compare_line:
			break

		#set up dictionary for the score terms from the compare file
		compare_conf_dict = {}
		term_counter = 0
		for item in compare_line.split(","):
			score_term_names.append(item.strip())
			#add the term to the dictionary
			compare_conf_dict[score_term_names[term_counter]] = item.strip()

			term_counter = term_counter + 1

		#get the delta between like terms (besides file) for the initial and compare, and determine if the compare is an improvement
		term_counter = 0
		improved_term_counter = 0

		#store all data in a list to be held on to in case this is one of the best improvements upon the original
		full_line_data = []

		for i in range(len(score_term_names_delta_write_file)):
			#write the comparisons
			#handle "file" differently since it is a string
			if "file" in score_term_names_delta_write_file[i]:
				#write the relevant file
				if "initial" in score_term_names_delta_write_file[i]:
					delta_data_file.write(initial_conf_dict["file"] + ",")
					full_line_data.append(initial_conf_dict["file"])
				if "compare" in score_term_names_delta_write_file[i]:
					delta_data_file.write(compare_conf_dict["file"] + ",")
					full_line_data.append(compare_conf_dict["file"])

				continue

			#likewise, handle "ligand" differently since it is a string (and may only exist in one file, if any)
			if "ligand" in score_term_names_delta_write_file[i]:
				#write the relevant ligand

				delta_data_file.write(initial_conf_dict["file"] + ",")
				full_line_data.append(initial_conf_dict["file"])

				continue


			#handle the other comparison terms

			#initial/compare/delta off the term for accessing in the corresponding dictionary
			stripped_term = score_term_names_delta_write_file[i].split("_initial")[0].split("_compare")[0].split("_delta")[0]

			if "initial" in score_term_names_delta_write_file[i]:
				delta_data_file.write(str(initial_conf_dict[stripped_term]) + ",")
				full_line_data.append(initial_conf_dict[stripped_term])
			if "compare" in score_term_names_delta_write_file[i]:
				delta_data_file.write(str(compare_conf_dict[stripped_term]) + ",")
				full_line_data.append(compare_conf_dict[stripped_term])
			#if delta
			if "delta" in score_term_names_delta_write_file[i]:

				#delta is initial minus compare
				delta_value = float(initial_conf_dict[stripped_term]) - float(compare_conf_dict[stripped_term])



				delta_data_file.write(str(delta_value) + ",")
				full_line_data.append(delta_value)




				#ddg,total_motifs,significant_motifs,real_motif_ratio,hbond_motif_count,hbond_motif_energy_sum,closest_autodock_recovery_rmsd,closest_autodock_recovery_ddg,strain_energy,total

				#delta is initial - compare
				#negative delta is improvements for: total_motifs,significant_motifs,real_motif_ratio,hbond_motif_count,
				#positive delta is improvements for: ddg, hbond_motif_energy_sum,closest_autodock_recovery_rmsd,closest_autodock_recovery_ddg,strain_energy,total(?)

				#we are primarily concerned with ddg, total_motifs, significant_motifs,real_motif_ratio, and hbond_motif_count
				#since ddg is the only positve delta improvement, if we take the negative of the values, it aligns with negative
				#could probably turn this into an argument at some point if we ever wanted to focus on other terms
				terms_of_interest = ["ddg","total_motifs","significant_motifs","real_motif_ratio", "hbond_motif_count"]

				if stripped_term == "ddg":
					delta_value = delta_value * -1

				#if the stripped term is one of the terms we want to focus on, increment the term counter and determine if the compare placement is an improvement
				if stripped_term in terms_of_interest:
					#increment the term counter and determine whether the compare placement is an improvement
					term_counter = term_counter + 1

					if delta_value < 0:
						improved_term_counter = improved_term_counter + 1

			#improvement ratio
			if "improvement_ratio" == score_term_names_delta_write_file[i]:
				#write the improvement ratio
				improvement_ratio = improved_term_counter / term_counter

				delta_data_file.write(str(improvement_ratio) + ",")
				full_line_data.append(str(improvement_ratio))

			#placement rmsd calc
			if "rmsd" == score_term_names_delta_write_file[i]:
				#open the compare placement pdb file, accessed with compare_conf_dict["file"]
				compare_placement_file = open(compare_conf_dict["file"],"r")

				#read the file and get the atom data
				compare_conf_atom_coords = {}


				for line2 in compare_placement_file.readlines():
					#HETATM lines have the ligand data
					if line2.startswith("HETATM "):
						#3rd position (index value 2 by splitting)
						unique_atom = line2.split()[2]

						#continue if Hydrogen
						#easy solution is if the first character is H, since n oother relevant element starts with H
						if unique_atom[0] == "H":
							continue

						#otherwise add the atom to the dictionary with xyz values as list
						compare_conf_atom_coords[unique_atom] = [float(line2.split()[6]),float(line2.split()[7]),float(line2.split()[8])]

				compare_placement_file.close()

				#iterate over the initial molecule atoms dict by keys and get the rmsd of heavy atoms

				#hold the distance sum and number of atoms
				distance_sum = 0
				num_atoms = 0

				#handling if for some reason atoms don't line up
				mismatch = False

				for atom in initial_conf_atom_coords.keys():
					
					if atom not in compare_conf_atom_coords.keys():
						mismatch = True
						break

					num_atoms = num_atoms + 1
					#get the distance between the two atoms
					atom_atom_distance = ((initial_conf_atom_coords[atom][0]-compare_conf_atom_coords[atom][0])**2 + (initial_conf_atom_coords[atom][1]-compare_conf_atom_coords[atom][1])**2 + (initial_conf_atom_coords[atom][2]-compare_conf_atom_coords[atom][2])**2) ** 0.5

					#apply the distance to the distance_sum
					distance_sum = distance_sum + atom_atom_distance

				#calculate the rmsd
				if num_atoms > 0:
					rmsd = distance_sum / num_atoms
				else:
					rmsd = 100

				if mismatch:
					rmsd = 100

				#write the rmsd
				delta_data_file.write(str(rmsd) + ",")
				full_line_data.append(str(rmsd))

		#write a newline
		delta_data_file.write("\n")

		#add the full line data to the all data
		all_placements.append(full_line_data)

		#determine if the improvement has a better improvement ratio than what is previously recorded
		if float(full_line_data[len(full_line_data) - 1]) > best_improvement_ratio:
			best_improvement_ratio = float(full_line_data[len(full_line_data) - 1])



	#look through all placements and write the ones with the best improvement ratio to the best file
	for placement in all_placements:
		if float(placement[len(placement) - 1]) >= best_improvement_ratio:
			#write the placement information to the best file
			for item in placement:
				delta_data_file_most_improved.write(str(item) + ",")

			#write newline
			delta_data_file_most_improved.write("\n")






