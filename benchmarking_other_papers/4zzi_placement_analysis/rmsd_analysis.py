#simple specific script that looks at the reference Gryniukova_Z26395438_placement_atoms_remapped.pdb (has atom names mapped to Rosetta palcements for RMSD comparison), and compares the reference to Rosetta palcements
#an output written with the ranked rmsd of placements

#imports
import os,sys

#decompress the lig_only_Z26395438_placements.tar.gz directory so we can use it
os.system("tar -xzf lig_only_Z26395438_placements.tar.gz")

#create a list that holds the rosetta placements and corresponding rmsd
rmsds_list = []

#iterate over the reference and get the atom coordinates, saved to a dictionary

ref_atoms = {}

#open a read stream
ref_in = open("Gryniukova_Z26395438_placement_atoms_remapped.pdb", "r")
for line in ref_in.readlines():
	atom_name = line.strip().split()[2]
	x = float(line.strip().split()[6])
	y = float(line.strip().split()[7])
	z = float(line.strip().split()[8])
	ref_atoms[atom_name] = [x,y,z]
ref_in.close()

#iterate over the rosetta placement files
for r,d,f in os.walk("lig_only_Z26395438_placements"):
	for file in f:
		if file.endswith(".pdb"):
			#set the total distance to 0
			total_distance = 0

			#set the atom counter to 0
			atom_counter = 0

			#get the molecule atom coordinates, and then compare to the reference
			test_in = open(r + "/" + file,"r")
			
			for line in test_in.readlines():
				#extract the atom name and coordinates
				atom_name = line.strip().split()[2]
				x = float(line.strip().split()[6])
				y = float(line.strip().split()[7])
				z = float(line.strip().split()[8])

				#derive the atom-atom distance
				a_a_distance = ((x - ref_atoms[atom_name][0])**2 + (y - ref_atoms[atom_name][1])**2 + (z - ref_atoms[atom_name][2])**2)**0.5

				total_distance = total_distance + a_a_distance
				atom_counter = atom_counter + 1

			#derive the rmsd as the total over the atom count
			rmsd = total_distance / atom_counter

			#add the rmsd to the rmsd list
			rmsds_list.append([file,rmsd])

sorted_list = sorted(rmsds_list, key=lambda x: x[1])

for item in sorted_list:
	print(item)

#end cleanup
os.system("rm -drf lig_only_Z26395438_placements")