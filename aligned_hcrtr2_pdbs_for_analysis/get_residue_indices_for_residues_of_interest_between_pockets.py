#the purpose of this script is to create a key for residues of interest between systems of interest (4s0v against 7l1u and 7sr8)

#import packages
import os,sys

#read in the list of residues of interest from 4s0v
#store to a list of tuples where index 0 is the residue 3 letter code and index 1 is the residue index
residues_of_interest = []

#read in the list
pocket_residue_file = open("4s0v_pocket_residues_unique.txt", "r")
for line in pocket_residue_file.readlines():
	#strip the newline
	stripped_line = line.strip()
	#make sure we are looking at a good line
	if len(stripped_line.split()) == 4:
		#then assign

		residue = stripped_line.split()[2]
		index = stripped_line.split()[3]
		residues_of_interest.append([residue,index])

#sanity check print of the residues of interest
#also make a list of the residue plus index for output

residues_of_interest_single_entry = []

for index in residues_of_interest:
	print(index)
	residues_of_interest_single_entry.append(str(index[0]) + "_" + str(index[1]))


#indexing should be the same for all atoms

#create dictionary that holds all atom data for each atom from residues of interest from each system
systems_atom_data = {}

#now, read through each local pdb file and get the rmsd of the residues of interest against 4s0v's residues
#the call to os.getcwd() means that this script is intended to be ran in ginsparg_thymelab_thesis/aligned_hcrtr2_pdbs_for_analysis
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		#make sure the file is in the working location and the file is a pdb file and that aligned is in the file name
		if r == os.getcwd() and file.endswith("_aligned.pdb"):

			#only work with 7sr8 and 7l1u
			if "7sr8" not in file and "7l1u" not in file and "4s0v" not in file:
				continue

			#extract the file name
			file_base_name = file.split("_")[0]
			print(file_base_name)

			#prepare systems dictionary for this system
			#add a dictionary using the base name as a key, and this nested dictionary will use the residue indices as keys
			systems_atom_data[file_base_name] = {}

			#read through the file and get the atom data from the residues of interest
			curr_file = open(file,"r")

			for line in curr_file.readlines():
				#only work with lines starting with ATOM
				if line.startswith("ATOM") == False:
					continue

				#check if the line corresponds to an residue index of interest
				#note, this method has potential to fail if the residue indices are 1000 or more, as chain id and residue index are no longer separated by a space
				#for this specific study, this easier method will work, but would need to be reworked if dealing with indices this high
				curr_index = line.split()[5]

				#print(line,curr_index)

				is_intersting = False
				#look through indices of interest
				for this_ind in residues_of_interest:
					if curr_index == this_ind[1]:
						#note is interesting and we will process
						is_intersting = True

				#do not process if not interesting
				if is_intersting == False:
					continue

				#skip if hydrogen (easy check is last character when splitting)
				stripped_line = line.strip()
				element = stripped_line.split()[len(stripped_line.split()) - 1]
				if element == "H":
					continue

				#print(line)

				#get the residue data
				#if this is the first atom, create a new dictionary entry
				if curr_index not in systems_atom_data[file_base_name].keys():
					systems_atom_data[file_base_name][curr_index] = {}

				#write the atom as a key (values should be unique) and coordinates as values
				atom_name = line.split()[2]
				atom_x = line.split()[6]
				atom_y = line.split()[7]
				atom_z = line.split()[8]

				systems_atom_data[file_base_name][curr_index][atom_name] = [float(atom_x),float(atom_y),float(atom_z)]

#system atom data should be filled now

#make a 2d dictionary that holds the rmsd between 4s0v and other systems that is the system on 1 dimension and the residue indices on the other
rmsd_dict = {}

#compare distances
for key in systems_atom_data.keys():
	#test print
	#print(systems_atom_data[key])
	#print(len(systems_atom_data))

	#skip 4s0v to not compare against self
	#if key == "4s0v":
	#	continue

	#add system dictionary to rmsd dictionary if it is not already present
	if key not in rmsd_dict.keys():
		rmsd_dict[key] = {}

	#iterate over each residue
	for index_key in systems_atom_data[key].keys():
		#variable to hold the distance sum
		distance_sum = 0

		#variable to hold the atom count (to get a mean to normalize for atom count)
		atom_count = 0

		#iterate over each atom for the given index for current system and 4s0v
		for atom_key in systems_atom_data[key][index_key].keys():
			#calculate the atom-atom distance
			#square root of sum of differences in x y z
			atom_atom_distance = ((systems_atom_data[key][index_key][atom_key][0] - systems_atom_data["4s0v"][index_key][atom_key][0])**2 + (systems_atom_data[key][index_key][atom_key][1] - systems_atom_data["4s0v"][index_key][atom_key][1])**2 + (systems_atom_data[key][index_key][atom_key][2] - systems_atom_data["4s0v"][index_key][atom_key][2])**2)**0.5

			#add the distance to the distance sum
			distance_sum = distance_sum + atom_atom_distance

			atom_count = atom_count + 1

		#with the distance sum, incorporate it into the rmsd dictionary
		rmsd_dict[key][index_key] = index_key
		print(key,index_key)

#write the rmsd dictionary to a csv file

#open write file
out_file = open("4s0v_pocket_residue_distances.csv", "w")

#lead with empty entry so first column can be indices
out_file.write(",")

#write header line of systems to first line
#hardcoding to just use 4s0v
for key in rmsd_dict["4s0v"].keys():
	
	for index in residues_of_interest_single_entry:
		split_index = index.split("_")[1]
		#write the matching index that also notes the residue
		if key == split_index:
			out_file.write(index + ",")

#cap first line
out_file.write("\n")

#iterate over keys
for key in rmsd_dict.keys():

	#variable to determine if the index has been written already
	#index_written = False

	out_file.write(key + ",")

	for index_key in rmsd_dict[key].keys():
		"""
		#first try to write index if it has not been written yet
		if index_written == False:
			index_written = True
			#pull the index from residues_of_interest_single_entry

		"""

		#now, write the corresponding rmsd for the residue at the index for the system
		out_file.write(str(rmsd_dict[key][index_key]) + ",")

	#cap line with newline
	out_file.write("\n")