import os,sys

#small script for making pymol session of all agonist results
#generally not meant to be reused, as other people's organizational structures can differ

#highest directory number (this also assumes that there are directories from 0 to this number)
max_dir_num = int(sys.argv[1])

#working directory
working_dir = sys.argv[2]

#ligand residue index
lig_index = str(sys.argv[3])

#string for residues to highlight
#example: 227+86+253+63+257
highlight_res_string = str(sys.argv[4])

#residue key file
residue_key_file = sys.argv[5]

for i in range(0,max_dir_num + 1):
	#move to directory
	os.chdir(str(i))

	#run the make session script
	os.system("python  /data/user/abgvg9/ginsparg_thymelab_thesis/discovery_placement_filtering/placement_pymol_session_creation/make_pymol_sessions_of_placements_pymol2.py " + working_dir + "/" + str(i) + " " + lig_index + " " + highlight_res_string + " " + residue_key_file)

	#make the folder to hold all sessions
	os.system("mkdir " + working_dir + "/placement_pymol_sessions")

	#copy the created sessions to the central location
	os.system("cp hbonds_only_session.pse " + working_dir + "/placement_pymol_sessions/" + str(i) + "_hbonds_only_session.pse")
	os.system("cp all_proteins_session.pse " + working_dir + "../placement_pymol_sessions/" + str(i) + "_all_proteins_session.pse")

	#go up a level
	os.chdir("..")
