#The purpose of this script is to look at all sdf files in the directory that this script is called from
#the script will real through each sdf file and get the MW of the molecule, counting only heavy atoms (non-hydrogen)
#The script will output a csv file that contains the MW of each file as well as a the average MW of all ligands

#import os and sys
import os,sys

#dictionary of elements potentially found in ligands (safe with Se and Si) to access their molecular weights
elements_atomic_weights = {
    'O': 15.999, 'N': 14.007, 'P': 30.974, 'S': 32.06, 'C': 12.011,
    'Se': 78.971, 'Si': 28.085, 'F': 18.998, 'Cl': 35.45, 'Br': 79.904,
    'I': 126.90
}


#get the working location
location = os.getcwd()

#declare trackers for mw sum and count of ligands to derive average
num_ligands = 0
total_mw = 0

#create a file to write all mw to
mw_file = open("ligand_mw.csv", "w")

#look over the location and get the MW of all sdf files
for r,d,f in os.walk(location):
	for file in f:
		#look only at local sdf files
		if r == location and file.endswith(".sdf"):

			#at a new file
			print(file)

			#variable to hold the ligand mw of heavy atoms
			lig_mw = 0

			#read through the file and get the mw
			lig_file = open(file, "r")

			for line in lig_file.readlines():
				#check if the line has at least 4 entries when split by spaces
				if len(line.split() > 4):
					#check if the entry at line.split()[3] is in the dictionary
					element = line.split()[3]

					if element in elements_atomic_weights.keys():
						#we have an atom, add the value to the ligand mw
						lig_mw = lig_mw + elements_atomic_weights[element]

			#done reading the file and we have its mw
			print(lig_mw)

			#add the mw to the total and increment the number of ligands
			num_ligands = num_ligands + 1
			total_mw = total_mw + lig_mw

			#write the ligand mw to the csv file
			mw_file.write(file + "," + str(lig_mw) + "\n")

#done looking at all ligands

#get the average mw and write it to the csv
average_mw = total_mw / num_ligands

mw_file.write("AVERAGE," + str(average_mw) + "\n")