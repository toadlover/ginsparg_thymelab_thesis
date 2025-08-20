#simple script to prepare and run rosetta on Cheaha for a single system
#very specific usage that is just thrown together for a quick batch of discovery, but I'll keep the logic and code here
#this is also specifically designed for the 7sr8 system

#imports
import os,sys

#get the system folder (i.e. /data/project/thymelab/august_2025_64b_agonist_results_pulldown/umass_conformer_receiving_and_round2_discovery/discovery_directories/m_276452____27418004____14224922____21064164__0)
#the folder will have a compressed test_params directory
system = sys.argv[1]

#enter the directory
os.chdir(system)

#unzip the test_params directory
os.system("tar -xzf test_params.tar.gz")

#write an args file
arg_file = open("args", "w")
arg_file.write("-constant_seed 1\n")
arg_file.write("-s /home/abgvg9/ginsparg_thymelab_thesis/7sr8_investigation/7sr8_receptor_only.pdb\n")
arg_file.write("-params_directory_path test_params/\n")
arg_file.write("-motif_filename /home/abgvg9/ginsparg_thymelab_thesis/7sr8_investigation/7sr8_pro_tak925_motif_S_O_N.motifs\n")
arg_file.write("-protein_discovery_locus 82\n")
arg_file.write("-constrain_relax_to_start_coords\n")
arg_file.write("-fa_atr_cutoff = -2\n")
arg_file.write("-fa_rep_cutoff = 150\n")
arg_file.write("-ddg_cutoff = -9\n")
arg_file.write("-best_pdbs_to_keep = 0\n")
arg_file.write("-binding_pocket_center_sf 122,125,160\n")
arg_file.write("-binding_pocket_dimensions_sf 4,5,5\n")
arg_file.write("-space_fill_cutoff_differential_score_sub 0.10\n")
arg_file.write("-collect_motifs_from_placed_ligand true\n")
arg_file.write("-significant_residues_for_motifs 85,82,62\n")
arg_file.write("-minimum_motifs_formed_cutoff 6\n")
arg_file.write("-minimum_significant_motifs_formed_cutoff 1\n")
arg_file.write("-mandatory_residues_for_motifs 85,82\n")
arg_file.write("-check_if_ligand_motifs_match_real true\n")
arg_file.write("-duplicate_dist_cutoff 1.2\n")
arg_file.write("-duplicate_angle_cutoff 0.45\n")
arg_file.write("-post_highresdock_fa_atr_rep_score true\n")
arg_file.write("-highresdock_with_whole_score_fxn true\n")
arg_file.write("-output_motifs_as_pdb false\n")
arg_file.close()

#run rosetta
os.system("/data/project/thymelab/rosetta_copy_cleaning_for_pr/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @args")

