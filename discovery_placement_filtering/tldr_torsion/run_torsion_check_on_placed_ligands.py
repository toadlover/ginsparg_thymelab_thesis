#note, you need to activate a conda environment that contains openbabel and rdkit for this to work
#on my machine, current environments that have them are my-rdkit-env and cpp_env

import os,sys

#the path to the location of files
#this program will look down at all directories within this location and perform the torsion calculation on each pdb
files_location = sys.argv[1]

#if files_location does not end with a backslash, append one for consistency
if files_location.endswith("/") == False:
	files_location = files_location + "/"


#the path to the tldr STRAIN torsion script, since their script does not appear to behave unless you are at the location of the script due to file dependencies
tldr_location = sys.argv[2]

#if tldr_location does not end with a backslash, append one for consistency
if tldr_location.endswith("/") == False:
	tldr_location = tldr_location + "/"

#go to the location
os.chdir(files_location)

#do os.walk on desired location
for r,d,f in os.walk(files_location):
	#in the location, look at all pdb placements
	for file in f:

		#select for .pdb files
		if file.endswith(".pdb"):

			#go to the file's location
			os.chdir(r)

			#for the pdb placement, extract the ligand from the pdb via grep
			#extract the name from the file type (i.e. remove the .pdb)
			file_no_suffix = file.split(".")[0]

			#extract the HETATM to a new pdb of only the ligand
			os.system("grep HETATM " + file + " > " + file_no_suffix + "_lig.pdb")

			#convert the isolated ligand to mol2 format
			os.system("obabel " + file_no_suffix + "_lig.pdb -O " + file_no_suffix + "_lig.mol2")

			#change the directory to where the torsion script is ??? (I don't think the scripts likes not being run outside of locally)
			os.chdir(tldr_location)

			#run the torsion script on the single placement mol2
			os.system("python Torsion_Strain.py " + r + "/" + file_no_suffix + "_lig.mol2")

			#extract the derived torsion energies (high and low) and rank the ligand
			#cat the contents to a master mol2 file at files_location
			os.system("cat " + r + "/" + file_no_suffix + "_lig_Torsion_Strain.csv >> " + files_location + "all_torsion_data.csv")
