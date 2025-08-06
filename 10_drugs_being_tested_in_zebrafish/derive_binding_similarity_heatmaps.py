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
			os.system("rm -drf motifs")
			os.system("mkdir motifs")
			os.chdir("motifs")

			os.system("/pi/summer.thyme-umw/2024_intern_lab_space/rosetta/source/bin/identify_ligand_motifs.linuxgccrelease @" + r + "/" + dire + "/args")


