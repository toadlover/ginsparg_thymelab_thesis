#simple specific script that looks at the 10 reference ligands placed in 5c8s by Singh, Li et al. in the Structure-Based Discovery of Inhibitors of the SARS-CoV-2 Nsp14 N7-Methyltransferase publication
#an output written with the ranked rmsd of placements for each ligand

#run as: python rmsd_analysis.py

#imports
import os,sys

#iterate over each remapped ligand pdb file in 5c8s_publication_placements
for r,d,f in os.walk("5c8s_publication_placements"):
	for file in f:
		#if the file is a remapped ligand file, work on doing an rmsd comparison
		if file.endswith("_remapped.pdb"):

			#derive the ligand name for later file naming and accessing the right session
			lig_name = file.split("_remapped.pdb")[0]

			#create a list that holds the rosetta placements and corresponding rmsd
			rmsds_list = []

			#iterate over the reference and get the atom coordinates, saved to a dictionary
			ref_atoms = {}

			#open a read stream
			ref_in = open(r + "/" + file, "r")
			for line in ref_in.readlines():
				atom_name = line.strip().split()[2]
				x = float(line.strip().split()[6])
				y = float(line.strip().split()[7])
				z = float(line.strip().split()[8])
				ref_atoms[atom_name] = [x,y,z]
			ref_in.close()

			#unzip the corresponding tar
			os.chdir(lig_name)
			os.system("tar -xzf lig_only.tar.gz")
			os.chdir("..")

			#iterate over the rosetta placement files
			for r2,d2,f2 in os.walk(lig_name + "/lig_only"):
				for file2 in f2:
					if file2.endswith(".pdb"):
						#set the total distance to 0
						total_distance = 0

						#set the atom counter to 0
						atom_counter = 0

						#get the molecule atom coordinates, and then compare to the reference
						test_in = open(r2 + "/" + file2,"r")
						
						print(r2 + "/" + file2)

						#skip if not a placement pdb (since the original empty pdb is also present)
						if "_conf" not in file2:
							continue

						for line in test_in.readlines():
							print(line)

							#extract the atom name and coordinates
							atom_name = line.strip().split()[2]

							#if atom is hydrogen, skip
							if atom_name.startswith("H"):
								continue

							#if the atom is otherwise not in, print out the issue (it means I messed up remapping)
							if atom_name not in ref_atoms.keys():
								print(atom_name,file2)

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
						rmsds_list.append([file2,rmsd])

			sorted_list = sorted(rmsds_list, key=lambda x: x[1])

			#for item in sorted_list:
			#	print(item)

			#write the sorted list to a file
			write_file = open(lig_name + "_rosetta_placement_rmsds.csv","w")
			for item in sorted_list:
				print(item)
				write_file.write(str(item[0]) + "," + str(item[1]) + "\n")



			#end cleanup
			os.system("rm -drf " + lig_name + "/lig_only")