#the purpose of this script is to operate in the working directory and look at all ligand placement pdbs
#this is a very niche script that is probably only being used once, so the logic is specific, but could later be expanded
#there should be a placement file ending with 0.pdb, which has the ddg (this is the really specific part)
#you then look for a file with the same name (remove .pdb), with _trimmed_real_motif_data.pdb added to the end, which has the real motif ratio
#this is done because of nuance with discovery with a single motif, and so you can't get a real motif ration originally, since you only use a single motif for discovery
#the naming schema of the 64B library is also hard to work with due to inconsistent underscore usage in ligand names
#we are goign to go ahead and run tldr strain on these too

#imports
import os,sys

#open a file to write the pdb data to
data_file = open("placements_data.csv", "w")

#write header
data_file.write("file,ligand,ddg,total_interactions,real_motifs,real_motif_ratio,strain\n")

working_location = os.getcwd()

#look over the current directory and look for relevant initial files
for r,d,f in os.walk(working_location):
	for file in f:
		if file.endswith("0.pdb") and working_location == r:
			
			os.chdir(working_location)

			#get the full file name
			full_file = r + "/" + file

			file_prefix = file.split(".pdb")[0]

			data_file.write(full_file + ",")

			system_ddg = "1000"

			#extract the ddg
			read_file = open(file, "r")
			for line in read_file.readlines():
				if "Post-HighResDock system ddG" in line:
					system_ddg = line.strip().split()[len(line.strip().split()) - 1]

			data_file.write(str(system_ddg) + ",")

			#now, extract the interaction/motifs data
			interactions_file = file_prefix + "_trimmed_real_motif_data.pdb"


			interactions = "0"
			real_motifs = "0"
			real_motif_ratio = "0"

			read_file = open(interactions_file, "r")
			for line in read_file.readlines():
				if "Real motif count" in line:
					real_motifs = line.strip().split()[len(line.strip().split()) - 1]
				if "Real motif ratio" in line:
					real_motif_ratio = line.strip().split()[len(line.strip().split()) - 1]
				if "Total motifs made" in line:
					interactions = line.strip().split()[len(line.strip().split()) - 1]

			data_file.write(str(interactions) + ",")
			data_file.write(str(real_motifs) + ",")
			data_file.write(str(real_motif_ratio) + ",")

			#finally get strain
			#corresponding mol2 files has the base followed by _lig.mol2
			#add the root to the file too (should be in same location)
			lig_file = r + "/" + file_prefix + "_lig.mol2"

			#move to the location of tldr strain (also need rdkit installed)
			os.chdir("/data/user/abgvg9/STRAIN/STRAIN_FILTER")

			#run the strain script on the mol2
			os.system("python Torsion_Strain.py " + lig_file)

			#this should make a torsion score csv file in the working directory
			#return to working directory
			os.chdir(working_location)

			#read the generated csv file to get the torsion strain upper and lower bounds (and take the average to use)
			torsion_data_file = r + "/" + file_prefix + "_lig_Torsion_Strain.csv"

			read_file = open(torsion_data_file,"r")

			torsion_average = "10"

			#should only be a single line
			for line in read_file.readlines():
				torsion_max = float(line.split(",")[1])
				torsion_min = float(line.split(",")[2])
				torsion_average = (torsion_max + torsion_min)/2

			data_file.write(str(torsion_average) + "\n")




