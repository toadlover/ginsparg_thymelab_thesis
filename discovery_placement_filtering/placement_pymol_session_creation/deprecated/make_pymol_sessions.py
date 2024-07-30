import os,sys

#small script for making pymol session of all agonist results
#generally not meant to be reused, as other people's organizational structures can differ

for i in range(0,37):
	#move to directory
	os.chdir(str(i))

	#run the make session script
	os.system("python  /data/user/abgvg9/ginsparg_thymelab_thesis/discovery_placement_filtering/placement_pymol_session_creation/make_pymol_sessions_of_placements_pymol2.py /scratch/abgvg9/discovery_results/top_1000_placement/agonist_12M_passing_placements/" + str(i) + " 282 227+86+253+63+257")

	#make the folder to hold all sessions
	os.system("mkdir ../placement_pymol_sessions")

	#copy the created sessions to the central location
	os.system("cp hbonds_only_session.pse ../placement_pymol_sessions/" + str(i) + "_hbonds_only_session.pse")
	os.system("cp all_proteins_session.pse ../placement_pymol_sessions/" + str(i) + "_all_proteins_session.pse")

	#go up a level
	os.chdir("..")
