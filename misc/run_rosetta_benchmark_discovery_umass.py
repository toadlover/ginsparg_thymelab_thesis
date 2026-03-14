#the purpose of this script is to run the rosetta benchmark on the dude library on the umass hpc

import os,sys

disc_out = "/home/ari.ginsparg-umw/patel_ginsparg_rosetta_lds_dude_benchmark/runs/3_13_26/disc_out"

os.chdir(disc_out)

#iterate over each folder, which is a system with an individual locus
for r,d,f in os.walk(disc_out):
	for dire in d:
		if r == disc_out:
			print(dire)

			os.chdir(disc_out)

			#break up the directory name to get the ssytem name and locus
			sys_name = dire.split("_")[0]
			locus = dire.split("_")[2]

			#move into the directory
			os.chdir(dire)

			#housekeeping, remove existing pdbs and other lingering files 
			os.system("rm -drf *pdb RO* *motifs *lig* *test_params* *txt")

			#copy the important files to the space (and just not have to deal with binding in the container execution)
			#pdb
			os.system("cp ../../../../files/all/" + sys_name + "/" + sys_name +  ".pdb . ")
			#motifs file
			os.system("cp ../../../../files/motifs_files/" + sys_name + "_ignorechain_FINAL_motifs_list_filtered_2_3_2023.motifs . ")
			#os.system("cp ../../../../files/FINAL_motifs_list_filtered_2_3_2023.motifs . ")
			#test_params
			os.system("cp -drf ../../../../files/all/" + sys_name + "/test_params . ")

			#write up an args file
			arg_file = open("args", "w")
			arg_file.write("-clean_pdb_name false\n")
			arg_file.write("-s " + sys_name + ".pdb\n")
			arg_file.write("-constant_seed 1\n")
			arg_file.write("-highresdock_with_whole_score_fxn false\n")
			arg_file.write("-fa_atr_rep_cutoff 10000000\n")
			arg_file.write("-fa_rep_cutoff 100\n")
			arg_file.write("-fa_atr_cutoff -5\n")
			arg_file.write("-ddg_cutoff 0\n")
			#arg_file.write("-motif_filename " + sys_name + "_ignorechain_FINAL_motifs_list_filtered_2_3_2023.motifs\n")
			arg_file.write("-motif_filename FINAL_motifs_list_filtered_2_3_2023.motifs\n")
			arg_file.write("-protein_discovery_locus " + locus + "\n")
			arg_file.write("-params_directory_path test_params/\n")
			arg_file.write("\n")
			arg_file.close()

			#run rosetta in a bsub job
			os.system("bsub -q short -W 8:00 -u \"\" -R \"rusage[mem=10000]\" \"singularity exec /pi/summer.thyme-umw/enamine-REAL-2.6billion/rosetta_condensed_6_25_2024.sif /rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @args && rm " + sys_name + ".pdb\"")
