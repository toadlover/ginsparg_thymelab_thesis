#script that runs autodock vina on placed ligands in a given directory
#goal is to determine if autodock is able to recover the predicted placement
#recovery will be determined via calculation of rmsd of the original rosetta placement versus all autodock placements

#need to have openbabel in your environment

import os,sys

#the path to the location of files
#this program will look down at all directories within this location and perform the torsion calculation on each pdb
files_location = sys.argv[1]

#if files_location does not end with a backslash, append one for consistency
if files_location.endswith("/") == False:
	files_location = files_location + "/"


#the path to the autodock executable
#include the entire path including the executable
#i.e. /data/user/abgvg9/autodock_vina_1_1_2_linux_x86/bin/vina
autodock_location = sys.argv[2]

#if tldr_location does not end with a backslash, append one for consistency
#if autodock_location.endswith("/") == False:
#	autodock_location = autodock_location + "/"

#go to the location
os.chdir(files_location)

#open output csv file to hold all data on autodock placements
data_csv = open("autodock_data.csv", "w")
#write header line
data_csv.write("system,model,rmsd,close_to_rosetta,vina_energy\n")

#do os.walk on desired location
for r,d,f in os.walk(files_location):
	#in the location, look at all pdb placements
	for file in f:
		#select for .pdb files
		if file.endswith(".pdb"):


			#go to the file's location
			os.chdir(r)

			#extract the ligand and convert it to .pdbqt

			#for the pdb placement, extract the ligand from the pdb via grep
			#extract the name from the file type (i.e. remove the .pdb)
			file_no_suffix = file.split(".")[0]

			#extract the HETATM to a new pdb of only the ligand
			os.system("grep HETATM " + file + " > " + file_no_suffix + "_lig.pdb")

			#convert the isolated ligand to pdbqt format
			os.system("obabel " + file_no_suffix + "_lig.pdb -O " + file_no_suffix + "_lig.pdbqt")

			#extract the receptor to a new pdb of only the receptor
			os.system("grep ATOM " + file + " > " + file_no_suffix + "_receptor.pdb")

			#convert the isolated receptor to pdbqt format, make it rigid using -xr flag
			os.system("obabel " + file_no_suffix + "_receptor.pdb -O " + file_no_suffix + "_receptor.pdbqt -xr")

			#derive the coordinates of all atoms in the ligand
			#used in rmsd calculation and for determining the centroid of the ligand to define the area for autodock to try to place
			#create a dictionary of non-hydrogen atoms with the atom with identifier as the key and a tuple of the x,y,z coordinates as the value
			ligand_atom_dict = {}

			#read the ligand pdbqt file to get the atom data
			original_lig_file = open(file_no_suffix + "_lig.pdbqt", "r")
			for line in original_lig_file.readlines():
				#if the line starts with "ATOM", it is an atom that we want to try to extract data
				if line.startswith("ATOM"):
					#make sure that the first character of split[2] isn't H (indicative of hydrogen; we shouldn't see any other chemical symbol that starts with H in organic chemistry)
					if line.split()[2][0] != "H":
						#set the atom entry in the dictionary
						ligand_atom_dict[line.split()[2]] = [float(line.split()[6]),float(line.split()[7]),float(line.split()[8])]

						#print(float(line.split()[6]),float(line.split()[7]),float(line.split()[8]))
			original_lig_file.close()

			#print(ligand_atom_dict)

			#derive the centroid and maximum dimensions of the ligand to be used in the autodock arg file

			#extrema values, defaulted to blank strings to be recognized as blank the first time they are encountered
			xmin = ""
			xmax = ""
			ymin = ""
			ymax = ""
			zmin = ""
			zmax = ""

			#iterate over each atom and define the extrema
			for atom in ligand_atom_dict.keys():
				
				#initial handling for fresh values
				if xmin == "":
					xmin = ligand_atom_dict[atom][0]
				#once it has a value, check to potentially set
				if xmin > ligand_atom_dict[atom][0]:
					xmin = ligand_atom_dict[atom][0]
								#initial handling for fresh values
				if xmax == "":
					xmax = ligand_atom_dict[atom][0]
				#once it has a value, check to potentially set
				if xmax < ligand_atom_dict[atom][0]:
					xmax = ligand_atom_dict[atom][0]

				#initial handling for fresh values
				if ymin == "":
					ymin = ligand_atom_dict[atom][1]
				#once it has a value, check to potentially set
				if ymin > ligand_atom_dict[atom][1]:
					ymin = ligand_atom_dict[atom][1]
								#initial handling for fresh values
				if ymax == "":
					ymax = ligand_atom_dict[atom][1]
				#once it has a value, check to potentially set
				if ymax < ligand_atom_dict[atom][1]:
					ymax = ligand_atom_dict[atom][1]

				#initial handling for fresh values
				if zmin == "":
					zmin = ligand_atom_dict[atom][2]
				#once it has a value, check to potentially set
				if zmin > ligand_atom_dict[atom][2]:
					zmin = ligand_atom_dict[atom][2]
								#initial handling for fresh values
				if zmax == "":
					zmax = ligand_atom_dict[atom][2]
				#once it has a value, check to potentially set
				if zmax < ligand_atom_dict[atom][2]:
					zmax = ligand_atom_dict[atom][2]

			#print(xmin, xmax, ymin, ymax, zmin, zmax)

			#define the centroid and dimensions for autodock to explore in about the centroid
			#centroid coordinates defined as 

			centroid = [xmax - ((xmax - xmin) / 2),ymax - ((ymax - ymin) / 2),zmax - ((zmax - zmin) / 2)]

			#define the distance to investigate in each dimension
			#defined as the difference of extrema multiplied by 1.5 to expand the area that autodock can investigate
			volume_lengths = [(xmax - xmin) * 1.5, (ymax - ymin) * 1.5, (zmax - zmin) * 1.5]

			#print(centroid)
			#print(volume_lengths)

			#write the arg file
			#arg_file = open(file_no_suffix + "_args.txt", "w")
			#arg_file.write("--receptor = " + r + "/" + file_no_suffix + "_receptor.pdbqt \n")
			#arg_file.write("--ligand = " + r + "/"  + file_no_suffix + "_lig.pdbqt \n")
			#arg_file.write("--center_x = " + str(centroid[0]) + "\n")
			#arg_file.write("--center_y = " + str(centroid[1]) + "\n")
			#arg_file.write("--center_z = " + str(centroid[2]) + "\n")
			#arg_file.write("--size_x = " + str(volume_lengths[0]) + "\n")
			#arg_file.write("--size_y = " + str(volume_lengths[1]) + "\n")
			#arg_file.write("--size_z = " + str(volume_lengths[2]) + "\n")
			#arg_file.write("--out = " + r + "/"  + file_no_suffix + "_autodock_placements.pdbqt \n")
			#arg_file.write("\n")
			#arg_file.close()

			#run autodock on the ligand
			#os.system(autodock_location + " " + file_no_suffix + "_args.txt")
			#print(autodock_location + " --receptor " + r + "/" + file_no_suffix + "_receptor.pdbqt --ligand " + r + "/"  + file_no_suffix + "_lig.pdbqt --center_x " + str(centroid[0]) + " --center_y " + str(centroid[1]) + " --center_z " + str(centroid[2]) + " --size_x " + str(volume_lengths[0]) + " --size_y " + str(volume_lengths[1]) + " --size_z " + str(volume_lengths[2]) + " --out " + r + "/"  + file_no_suffix + "_autodock_placements.pdbqt --log " + r + "/"  + file_no_suffix + "_autodock_log.txt --num_modes 25 --seed 0")
			os.system(autodock_location + " --receptor " + r + "/" + file_no_suffix + "_receptor.pdbqt --ligand " + r + "/"  + file_no_suffix + "_lig.pdbqt --center_x " + str(centroid[0]) + " --center_y " + str(centroid[1]) + " --center_z " + str(centroid[2]) + " --size_x " + str(volume_lengths[0]) + " --size_y " + str(volume_lengths[1]) + " --size_z " + str(volume_lengths[2]) + " --out " + r + "/"  + file_no_suffix + "_autodock_placements.pdbqt --log " + r + "/"  + file_no_suffix + "_autodock_log.txt --num_modes 25 --seed 0")

			#determine the rmsd of all placements against the rosetta placement and if any are close (within 2 angstroms)
			#read the placements pdbqt file
			#need to try this in case autodock gets no placements, which would crash this script
			try:
				placements_file = open(r + "/"  + file_no_suffix + "_autodock_placements.pdbqt", "r")
			except:
				#behavior if no placements were obtained
				#write a single line to the csv file that states that this ligand had no placements from autodock
				data_csv.write(file_no_suffix + ",0,10000,False,10000\n")
				continue

			#hold the model number
			model_num = ""

			#dictionary to hold the model numbers and the corresponding rmsd off of the rosetta palcements, autodock energy, and if they are within the expected rmsd (<2 angstroms)
			model_rmsd = {}

			#working dictionary to hold the atom coordinates of the placement
			model_ligand_atom_dict = {}

			#autodock energy score
			vina_energy = ""

			for line in placements_file.readlines():
				#get model number at start of model
				if line.startswith("MODEL"):
					model_num = line.split()[1].strip()

				#get atom data
				if line.startswith("ATOM"):
					#make sure that the first character of split[2] isn't H (indicative of hydrogen; we shouldn't see any other chemical symbol that starts with H in organic chemistry)
					if line.split()[2][0] != "H":
						#set the atom entry in the dictionary
						model_ligand_atom_dict[line.split()[2]] = [float(line.split()[6]),float(line.split()[7]),float(line.split()[8])]

				#get energy from autodock
				if line.startswith("REMARK VINA RESULT"):
					vina_energy = line.split()[3]

				#end of model behavior to calculate RMSD
				if line.startswith("ENDMDL"):
					#hold sum of atom distances
					atom_distance_sum = 0

					#hold sum of atom counts
					atom_count = 0

					#run through each atom to get the distance
					for atom in model_ligand_atom_dict.keys():
						#only attempt if atom is in both
						if atom in ligand_atom_dict.keys():
							#3d distance formula for the atom
							distance = (((model_ligand_atom_dict[atom][0] - ligand_atom_dict[atom][0]) ** 2) + ((model_ligand_atom_dict[atom][1] - ligand_atom_dict[atom][1]) ** 2) + ((model_ligand_atom_dict[atom][2] - ligand_atom_dict[atom][2]) ** 2)) ** 0.5

							#add the distance to the distance sum and increment the atom count
							atom_distance_sum = atom_distance_sum + distance
							atom_count = atom_count + 1

					#rmsd calculation
					rmsd = atom_distance_sum / atom_count

					close_placement = (rmsd <= 2)

					#add the data to the model_rmsd dictionary
					model_rmsd[model_num] = [str(rmsd), str(close_placement), vina_energy]

					#wipe the model ligand atom dictionary
					model_ligand_atom_dict = {}

			for placement in model_rmsd.keys():
				#print(model_rmsd[placement])
				#write the placement info to csv
				#data_csv.write("system,model,rmsd,close_to_rosetta,vina_energy\n")
				data_csv.write(file_no_suffix + "," + placement + "," + model_rmsd[placement][0] + "," + model_rmsd[placement][1] + "," + model_rmsd[placement][2] + "\n")




			#write the closest recovery rmsd to placement csv

			#cleanup
			#rm *lig.pdb* *receptor.pdb* *args.txt *log*txt *placements.pdb*
			#clean intermediates
			os.system("rm *lig.pdb* *receptor.pdb*")
