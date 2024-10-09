#the purpose of this script is to automate the usage of run_torsion_check_on_placed_ligands.py on all directories at and below whe input directory for this script
#this script takes in 3 arguments: argument 1 is for the location to operate from, argument 2 is for the path and script of run_torsion_check_on_placed_ligands.py, and argument 3 is for the paty to the tldr STRAIN script (just path, no script)
#i.e.: python run_torsion_on_directories_of_files.py /scratch/abgvg9/discovery_results/antagonists/placements_for_refinement/expanded_conformer_set_placements/low_ddg_best_refined_placements/ /scratch/abgvg9/STRAIN/STRAIN_FILTER

#this script requires an environment that has openbabel (obabel) and rdkit, since the torsion script calls for it     

import os,sys

#get the arguments
starting_directory = sys.argv[1]
if starting_directory.endswith("/") == False:
	starting_directory = starting_directory + "/"

torsion_script = sys.argv[2]
if "run_torsion_check_on_placed_ligands.py" not in torsion_script:
	print("You are not including the run_torsion_check_on_placed_ligands.py script, exiting now.")
	exit()

strain_directory = sys.argv[3]
if starting_directory.endswith("/") == False:
	starting_directory = starting_directory + "/"

#list to hold all directories that have been processed already
processed_directories = []

#move to starting location
os.chdir(starting_directory)

#walk over the starting directory, and iterate over all files 
for r,d,f in os.walk(starting_directory):
	for file in f:
		#if the file root is in processed_directories, then continue because we do not want to process the directory again
		if r in processed_directories:
			continue
		#make sure that this is a "compare" refined pdb file
		if "compare" in file and file.endswith(".pdb"):
			#add the root to processed_directories so that we don't process the directory that this file came from again
			processed_directories.append(r)
			print(r)

			#move to the root
			os.chdir(r)

			#make a folder called "strain" and copy all compare pdb files into the strain folder
			os.system("mkdir strain")
			os.system("cp *compare*.pdb strain")

			#move into the new directory
			os.chdir("strain")

			#run the torsion script
			os.system("python " + torsion_script + " " + r + " " + strain_directory)


"""
cd 5
mkdir strain
cp *compare*.pdb strain
cd strain
python /data/user/abgvg9/ginsparg_thymelab_thesis/discovery_placement_filtering/tldr_torsion/run_torsion_check_on_placed_ligands.py /scratch/abgvg9/discovery_results/antagonists/placements_for_refinement/expanded_conformer_set_placements/low_ddg_best_refined_placements/5/strain /scratch/abgvg9/STRAIN/STRAIN_FILTER
cd ../..
"""