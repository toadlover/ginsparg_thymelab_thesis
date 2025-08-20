#this script is intended to run Rosetta's determine_real_motifs protocol on a set of PDBs
#run this script on the directory that you are currently in

#imports
import os,sys

#run over the current directory and look at all pdb files, and runthe motifs collection
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		if file.endswith(".pdb") and r == os.getcwd():
			#work with the pdb

			file_base = file.split(".pdb")[0]

			#write an args file
			arg_file = open(file_base + "_args", "w")



			arg_file.write("-constant_seed 1\n")
			arg_file.write("-s " + file +" \n")
			arg_file.write("-motif_filename /data/user/abgvg9/FINAL_motifs_list_filtered_2_3_2023.motifs\n")
			arg_file.write("-collect_motifs_from_placed_ligand true\n")
			arg_file.write("-minimum_motifs_formed_cutoff 6\n")
			arg_file.write("-minimum_significant_motifs_formed_cutoff 1\n")
			arg_file.write("-check_if_ligand_motifs_match_real true\n")
			arg_file.write("-duplicate_dist_cutoff 1.2\n")
			arg_file.write("-duplicate_angle_cutoff 0.45\n")
			arg_file.write("-output_motifs_as_pdb false\n")
			arg_file.write("-in::ignore_unrecognized_res true\n")

			arg_file.close()

			#run the app on the file
			os.system("/data/project/thymelab/rosetta_copy_cleaning_for_pr/rosetta/source/bin/determine_real_motifs.linuxgccrelease @" + file_base + "_args")






