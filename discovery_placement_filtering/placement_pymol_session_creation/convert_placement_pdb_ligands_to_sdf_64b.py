#this script uses obabel to convert ligands in the working directory from a rosetta placement pdb file to sdf format
#modified version that handles ligand names for ligands from the 64b library

import os,sys

location = os.getcwd()

for r,d,f in os.walk(location):
	for file in f:
		#if the file is a pdb file and in the working location only
		if file.endswith(".pdb") and r == location:
			
			#print the file
			print(file)

			#pull the HETATM data into a new file based on the placement info (the 3rd to last, second to last, and last sections when splitting the file name by underscores)
			#lig_pdb = file.split("_")[len(file.split("_")) - 3] + "_" + file.split("_")[len(file.split("_")) - 2] + "_" + file.split("_")[len(file.split("_")) - 1]
			lig_pdb = file.split("only_")[1]

			#write the HETATM lines to this file
			os.system("grep HETATM " + file + " > " + lig_pdb)

			#make lig_sdf string
			lig_sdf = lig_pdb.split(".")[0] + ".sdf"

			#use openbabel to convert the ligand to a sdf
			os.system("obabel " + lig_pdb + " -O " + lig_sdf)

			#append the new sdf to an sdf of all ligands as well
			os.system("cat " + lig_sdf + " >> all.sdf")
