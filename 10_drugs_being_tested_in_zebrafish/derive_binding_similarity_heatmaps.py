#this is a script with a specific single-use pipeline purpose
#the purpose of this script is to go through the 10_drugs_being_tested_in_zebrafish folder and look at every pdb at the top level (or folder at this level if the pdbs have already been organized)
#each system needs to be aligned against each other for downstream analyses (and also assumes the same receptor being used, HCRTR2 in this case)
#for each system, Rosetta's identify_ligand_motifs script will be run on the system to collect motifs off the placed ligand
#from here, the motif PDBs that are generated will be used and mapped against a single system to note which residue was used across all systems (since different indexing exists)
#the AllMattMotifs.motifs file will be used to take the sum of packing and hbond energies at each residue (the value will be 0 for any residues found in other systems but not the working system)
#finally, a 2D grid will be made to compare the energy sums per residue in each system
#a second 2D grid will be made in the same format, but with binary 1/0 values to derive jaccard similarity of binding mode
#the grid can be used to make a heatmap, along with determining the jaccard similarity

#this script should be called from ginsparg_thymelab_thesis/10_drugs_being_tested_in_zebrafish, where it is located
#currently the call to identify_ligand_motifs is hard-coded for its location, but this can be changed if it needs to be run elsewhere (no plan to make a flag for now, since this is very specific)

#imports
import os,sys

#starting location
starting_location = os.getcwd()

#iterate over the current location to determine any pdb files in the current location
#for each pdb, make a new folder and move the pdb into the folder
#it is ok if an analysis has already been run and previous files were moved into created folders
#these folders can be reused for a rerun
for r,d,f in os.walk(starting_location):
	for file in f:
		#only look for pdb files at the top of the directory
		if r == starting_location and file.endswith(".pdb"):
			#make a directory and move the file into the directory

			file_base = file.split(".pdb")[0]

			os.system("mkdir " + file_base)
			os.system("mv " + file + " " + file_base)

#run over the 4s0v_aligned.pdb to get the center of mass of each residue. this will be used to map residues from motifs from other systems 
#make a list of 2 entry lists, with the first entry being the residue and the second being the coordinates of its center of mass (not including hydrogens)
reference_residues_com_list = []

#open the reference file
#again, this is hard-coded for HCRTR2 with a pre-aligned 4s0v
reference_file = open("../aligned_hcrtr2_pdbs_for_analysis/4s0v_aligned.pdb", "r")

#read the reference file
#hold the working residue
working_residue = -1
working_residue_name = ""
working_residue_atom_coords = []

#hold the list of atoms for the working residue

for line in reference_file.readlines():
	if line.startswith("ATOM"):
		residue_name = line[17:20].strip()
		residue_number = line[22:26].strip()

		#write over the default
		if working_residue == -1 and working_residue_name == "":
			working_residue = residue_number
			working_residue_name = residue_name

		#determine if we are on a new residue, and if so, we need to conclude the previous residue
		#conclude if there is a mismatch
		if residue_number != working_residue:
			#hold the sum of x,y,z coordinates and atom count, and derive the average of each for the residue center of mass
			x_sum = 0
			y_sum = 0
			z_sum = 0
			n_atoms = len(working_residue_atom_coords)

			for atom in working_residue_atom_coords:
				x_sum = x_sum + atom[0]
				y_sum = y_sum + atom[1]
				z_sum = z_sum + atom[2]

			#compose the list to add to the reference_residues_com_list
			residue_list = [residue_name + residue_number, [x_sum/n_atoms,y_sum/n_atoms,z_sum/n_atoms]]
			reference_residues_com_list.append(residue_list)

			print(residue_name + residue_number, [x_sum/n_atoms,y_sum/n_atoms,z_sum/n_atoms])

			#wipe the working residue list
			working_residue_atom_coords = []

			#assign the working residue number and name to the current
			working_residue = residue_number
			working_residue_name = residue_name


		x = float(line[30:38].strip())
		y = float(line[38:46].strip())
		z = float(line[46:54].strip())

		#add the coordinates as a list to working_residue_atom_coords
		working_residue_atom_coords.append([x,y,z])

#conclude the final residue after concluding the loop
#hold the sum of x,y,z coordinates and atom count, and derive the average of each for the residue center of mass
x_sum = 0
y_sum = 0
z_sum = 0
n_atoms = len(working_residue_atom_coords)

for atom in working_residue_atom_coords:
	x_sum = x_sum + atom[0]
	y_sum = y_sum + atom[1]
	z_sum = z_sum + atom[2]

#compose the list to add to the reference_residues_com_list
residue_list = [working_residue_name + working_residue, [x_sum/n_atoms,y_sum/n_atoms,z_sum/n_atoms]]
reference_residues_com_list.append(residue_list)
print(working_residue_name + working_residue, [x_sum/n_atoms,y_sum/n_atoms,z_sum/n_atoms])


#create a dictionary where each key is a system and the contents are a list with the mapped residue relative to 4s0v and the energy sum
system_residue_energies_dict = {}

#now, run through each folder at the starting level and attempt to run the motif collection and then map the residues against 4s0v
for r,d,f in os.walk(starting_location):
	for dire in d:
		#look at top level directories
		#check if there is a pdb file in the directory that matches. If so, we can work here
		if r == starting_location and os.path.exists(r + "/" + dire + "/" + dire + ".pdb"):
			print(r + "/" + dire + "/" + dire + ".pdb")



			#now, write an args file for the system and run Rosetta's identify_ligand_motifs
			arg_file = open(r + "/" + dire + "/args", "w")
			arg_file.write("-ignore_unrecognized_res\n")
			arg_file.write("-s " + r + "/" + dire + "/" + dire + ".pdb\n")
			arg_file.write("-hb_score_cutoff -0.3\n")
			arg_file.write("-water_score_cutoff -0.3\n")
			arg_file.write("-pack_score_cutoff -0.5\n")
			arg_file.close()

			#now, move into the directory, create another directory within to have the motifs pdbs and file be generated, and then run Rosetta
			#if there is an existing motifs folder, delete it
			os.chdir(dire)

			#check and sanitize the 3 letter code of the file if needed
			#Rosetta can not handle ligands with 3 letter codes that contain a -, so replace - with X
			#should be as simple as calling sed on the file
			os.system("sed -i 's/PV-/PVX/g' " + r + "/" + dire + "/" + dire + ".pdb")

			os.system("rm -drf motifs")
			os.system("mkdir motifs")
			os.chdir("motifs")

			os.system("/pi/summer.thyme-umw/2024_intern_lab_space/rosetta/source/bin/identify_ligand_motifs.linuxgccrelease @" + r + "/" + dire + "/args")

			#now, we can map the motif residues to the reference
			#iterate over each motif pdb that was derived, will be in the current location
			for r2,d2,f2 in os.walk(os.getcwd()):
				for motif_pdb in f2:
					if "_Ligatoms_" in motif_pdb and motif_pdb.endswith(".pdb"):
						print(motif_pdb)

						#derive the residue code and residue number
						motif_residue_code = motif_pdb[0:3]
						motif_residue_number = motif_pdb[3:6].strip("_")

						#prepare a list to hold atom coordinates so that the center of mass can be determined
						motif_residue_atom_coords = []

						#read the file to get the residue atom coordinate data
						read_motif_pdb = open(motif_pdb,"r")
						for line in read_motif_pdb:
							if line.startswith("ATOM"):
								#skip hydrogens
								if line.strip().endswith("H"):
									continue
								residue_name = line[17:20].strip()
								residue_number = line[22:26].strip()

								if residue_name != motif_residue_code:
									print("Warning, we somehow have a residue mismatch in the motif! " + motif_residue_code + " " + residue_name)
									continue

								x = float(line[30:38].strip())
								y = float(line[38:46].strip())
								z = float(line[46:54].strip())

								#add the coordinates as a list to working_residue_atom_coords
								motif_residue_atom_coords.append([x,y,z])

						#now that we have all atom coordinates, derive the center of mass
						x_sum = 0
						y_sum = 0
						z_sum = 0
						n_atoms = len(motif_residue_atom_coords)
						for atom in motif_residue_atom_coords:
							x_sum = x_sum + atom[0]
							y_sum = y_sum + atom[1]
							z_sum = z_sum + atom[2]

						#compose the list to add to the reference_residues_com_list
						residue_com = [x_sum/n_atoms,y_sum/n_atoms,z_sum/n_atoms]

						#now, iterate over all residues in the reference to determine which residue (of the same type) is the closest by COM distance
						#hold the closest residue as a list of the reference residue code+index and the distance
						closest_residue = ["XXXXXX",100000]

						for residue in reference_residues_com_list:
							#determine if the working residue is the same kind of residue, continue if not
							if residue[0].startswith(motif_residue_code):
								#derive the distance of COM
								distance = (((residue_com[0] - residue[1][0])**2)+((residue_com[1] - residue[1][1])**2)+((residue_com[2] - residue[1][2])**2))**0.5

								#if the distance is smaller than the current closest_residue distance, replace with current
								if distance < closest_residue[1]:
									closest_residue = [residue[0],distance]

						print(closest_residue)

						motif_energy_sum = 0

						#finally, extract the sum of packing and hbond scores from the motifs file
						motifs_file = open("AllMattMotifs.motifs","r")
						for line in motifs_file:
							#line must start with single
							if line.startswith("SINGLE"):
								#line must have the original (not mapped) residue name and code
								if (motif_residue_code + motif_residue_number) in line:
									#extract the packing and hydrogen bond scores
									packing_score = float(line.strip().split()[len(line.split()) - 1].split("Packing_score:")[1].split("_")[0])
									hbond_score = float(line.strip().split()[len(line.split()) - 1].split("Hbond_score:")[1])
									motif_energy_sum = packing_score + hbond_score

						#now, add the mapped residue with the motif energy sum to the system dictionary
						#if this is the first entry for this key, make a new list at the entry
						if dire not in system_residue_energies_dict.keys():
							system_residue_energies_dict[dire] = []

						#append the mapped residue and energy tuple to the list at the key
						system_residue_energies_dict[dire].append([closest_residue[0],motif_energy_sum])




			#at end, reset and move back up to the starting location
			os.chdir(starting_location)

			#now at end, test print of dictionary to confirm it looks good
			for key in system_residue_energies_dict.keys():
				print(system_residue_energies_dict[key])
