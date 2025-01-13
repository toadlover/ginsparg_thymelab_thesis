#the purpose of this script is to take the coordinates of select binding pocket residues and calculate the rmsd of the residue against 4s0v (after being previously aligned to 4s0v)

#this is intended to be a 1-off script to process data, but could be adapted in the future
#this script will produce a csv file that I can make a heatmap off of to determine the rmsd per residue per structure against 4s0v

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
#test print
for key in systems_atom_data.keys():
	print(systems_atom_data[key])
	print(len(systems_atom_data))



		