import os
import re
import os.path
import sys
import time

loc1 = str(sys.argv[1])
loc2 = str(sys.argv[2])
fa_rep_cutoff = str(sys.argv[3])
fa_atr_cutoff = str(sys.argv[4])
sub_area_score_cutoff = str(sys.argv[5])
print_space_fill_matrix = str(sys.argv[6])
sub_area_score_differential_cutoff = str(sys.argv[7])

os.mkdir("disc_out")
for f in os.listdir("/data/user/abgvg9/benchmark/motif_collection/files/neighbor_out"):
    d = f.replace("output_", "")
    d = d.replace(".txt", "")
    o = open("/data/user/abgvg9/benchmark/motif_collection/files/neighbor_out/" + f, "r")
    count = 0
    
    #added march 2024 by Ari, code to get a list of all residues in order to make a shape to represent a binding pocket for the space fill method (which did not exist when the benchmark was originally written)
    #open a second instance of the neighbor_out/output_*txt file to read all lines and get all positions
    out_file = open("/data/user/abgvg9/benchmark/motif_collection/files/neighbor_out/" + f, "r")
    #a string to build that is the indices (in string form)
    position_string = ""
    for line in out_file.readlines():
        #if "ROS" is in the line, we can extract a position
        #residue index should be index 3 in the line if split (make sure line split length is at least 3, ignor otherwise (it should always work))
        if "ROS" in line and len(line.split()) >= 3:
            #add the string of the index to the list, tail with a comma an no spaces
            #print(f,line,line.split()[3])
            position_string = position_string + str(line.split()[1]) +","

    #remove the tail comma
    position_string = position_string[:-1]

    #close the out_file stream
    out_file.close()

    #derive the center of the binding region and side length of area to investigate
    #derived from /benchmark/files/all/(ligname)/(ligname)-lig.pdb
    #derived as average of coordinates of extrema in xyz (centroid)
    #side length derived as 1/2 largest difference between min/max * 1.1
    lig_file = open("/data/user/abgvg9/benchmark/files/all/" + d + "/" + d + "-lig.pdb", "r")

    xmax = ""
    xmin = ""
    ymax = ""
    ymin = ""
    zmax = ""
    zmin = ""

    #read the lig file and find its centroid
    for line in lig_file.readlines():
        #ignore water lines
        if "HOH" in line:
            continue

        #x is line.split()[5], y is 6, z is 7
        if xmax == "" or int(line.split()[5].split(".")[0]) > xmax:
            xmax = int(line.split()[5].split(".")[0])

        if xmin == "" or int(line.split()[5].split(".")[0]) < xmin:
            xmin = int(line.split()[5].split(".")[0])

        if ymax == "" or int(line.split()[6].split(".")[0]) > ymax:
            ymax = int(line.split()[6].split(".")[0])

        if ymin == "" or int(line.split()[6].split(".")[0]) < ymin:
            ymin = int(line.split()[6].split(".")[0])

        if zmax == "" or int(line.split()[7].split(".")[0]) > zmax:
            zmax = int(line.split()[7].split(".")[0])

        if zmin == "" or int(line.split()[7].split(".")[0]) < zmin:
            zmin = int(line.split()[7].split(".")[0])

    #derive the centroid
    centroid = [int(xmax - ((xmax - xmin)/2)), int(ymax - ((ymax - ymin)/2)), int(zmax - ((zmax - zmin)/2))]

    #determine the longest difference
    side_length = xmax - xmin
    if ymax - ymin > side_length:
        side_length = ymax - ymin

    if zmax - zmin > side_length:
        side_length = zmax - zmin

    #multiply by 1.1
    #we may not want to do that
    side_length = side_length * 1.1
    
    #to conform to using binding_pocket_dimensions_sf flag, make list/vector for it
    x_length = str((xmax-xmin) / 2).split(".")[0]
    y_length = str((ymax-ymin) / 2).split(".")[0]
    z_length = str((zmax-zmin) / 2).split(".")[0]
    
    #convert to a string
    length_string = x_length + "," + y_length + "," + z_length

    #convert centroid and side length to strings
    centroid_str = str(centroid[0]) + "," + str(centroid[1]) + "," + str(centroid[2])
    side_length_str = str(side_length)

    print(d, centroid, centroid_str, length_string)

    m = {}
    for line in o:
        if count % 2 == 0:
            n = re.search('ROS (.*) Res', line).group(1)
            p = ""
            if n in m:
                os.mkdir("disc_out/" + d + "_res_" + n + "_" + str(m[n]))
                p = "disc_out/" + d + "_res_" + n + "_" + str(m[n])
                m[n] += 1
            else:
                os.mkdir("disc_out/" + d + "_res_" + n + "_0")
                p = "disc_out/" + d + "_res_" + n + "_0"
                m[n] = 1

            os.chdir(p)
            t = open(p.split("disc_out/", 1)[1] + "_args", "x")
            t.write(
"""-parser:protocol /data/user/abgvg9/benchmark/files/discovery_files/lite_enzdes.xml 

#CHANGE THIS
-s /data/user/abgvg9/benchmark/files/all/""" + d + """/""" + d + """.pdb

#used constant seed for the sake of consistency
-constant_seed 1

#optional flags to collect motifs off of placed ligands and see if the motifs resemble real ones
-collect_motifs_from_placed_ligand true
-check_if_ligand_motifs_match_real true
-duplicate_dist_cutoff 0.8
-duplicate_angle_cutoff 0.3

#do not keep motifs of placements to help save on space
-output_motifs_as_pdb false

#use the whole score function for ligand docking
-highresdock_with_whole_score_fxn false

#use whole ligand.wts score function filter
#-score_with_ligand_wts_function true
-fa_atr_rep_cutoff 10000000
#-ligand_wts_fxn_cutoff 10000

#3_22_2022, now using motif lists that remove homologous pdb motifs
#these files are now located in /data/user/abgvg9/benchmark/files named "d"_ignorechain_condensed_motif_list_08_06_21.motifs

#using 2022 motif list
-motif_filename """ + "/data/user/abgvg9/benchmark/files/" + d + "_ignorechain_FINAL_motifs_list_filtered_2_3_2023.motifs" """
#using 2021 motif list
#-motif_filename """ + "/data/user/abgvg9/benchmark/files/" + d + "_ignorechain_condensed_motif_list_08_06_21.motifs" """
#old
#-motif_filename """ + loc1 + """

#CHANGE THIS
-protein_discovery_locus """ + n + """

#space fill positions and score cutoff for sub-area
#-target_residues_sf  """ + position_string + """
-binding_pocket_center_sf """ + centroid_str + """
#-binding_pocket_radius_sf """ + side_length_str + """
-binding_pocket_dimensions_sf """ + length_string + """
-space_fill_cutoff_score_sub """ + sub_area_score_cutoff + """
-space_fill_cutoff_differential_score_sub """ + sub_area_score_differential_cutoff + """
#output the space fill matrix (probably won't always want this, as it slows things down and uses up extra memory)
-output_space_fill_matrix_pdbs """ + print_space_fill_matrix + """

#CHANGE THIS (path to test_params directory for receptor system)
-params_directory_path /data/user/abgvg9/benchmark/files/all/""" + d + """/test_params/

#-mute "Code tracing" core.init basic.random.init_random_generator core.chemical.GlobalResidueTypeSet protocols.motifs.motif_utils core.scoring.etable core.import_pose.import_pose core.conformation.Conformation core.pack.pack_missing_sidechains core.pack.task core.scoring.ScoreFunctionFactory basic.io.database core.pack.dunbrack.RotamerLibrary core.pack.pack_rotamers core.pack.interaction_graph.interaction_graph_factory core.io.pose_from_sfr.PoseFromSFRBuilder


-score::weights /data/user/abgvg9/benchmark/files/discovery_files/enzdes_2.5.wts
-dtest 7.5
-r2 1.1
-r1 4.5
-z1 0.1
-z2 1.0
-dump_motifs true
-output_file output.file
-motifs::data_file data.file
-motifs::minimize false
-minimize_dna false
-enzdes::run_ligand_motifs true
-patch_selectors SPECIAL_ROT
-special_rotweight -20.0
-rotlevel 8
-enzdes::detect_design_interface
-enzdes::cut1 4.0
-enzdes::cut2 6.0
-enzdes::cut3 6.0
-enzdes::cut4 6.0
-enzdes::design_min_cycles 1
-enzdes::chi_min
-in::ignore_unrecognized_res true
-enzdes::cstfile null.cst
-ex1
-ex2
-ex3
-ex4
-extrachi_cutoff 1
-out::file::o ./enz_score.out
-jd2:enzdes_out true
-soft_rep_design
-ex1aro
-ex1aro:level 6
-ex2aro
-ex2aro:level 6
-extrachi_cutoff 1
-soft_rep_design
-flip_HNQ
-nstruct 1
-enzdes::no_unconstrained_repack
-linmem_ig 10
-enzdes::lig_packer_weight 2.5
-num_repacks 4
-patch_selectors SPECIAL_ROT
-output_residue_energies
-run_motifs
-ignore_unrecognized_res
-ndruns 1
-fa_rep_cutoff """ + fa_rep_cutoff + """
-fa_atr_cutoff """ + fa_atr_cutoff + """
-out:prefix test_"""
            )
            t.close()
            t = open(p.split("disc_out/", 1)[1] + ".job", "x")
            t.write(
"""#!/bin/bash
#SBATCH -p medium # Partition to submit to
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1000 # Runtime in minutes
#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH -o hostname_%A_%a.out # Standard out goes to this file
#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n""" +

#need to change this line, depending on the name of the binary that needs to be called
#switching to using protocol, which will be integrated into main Rosetta code
"""/usr/bin/time -f "%e" """ + loc2 + """/source/bin/ligand_discovery_search_protocol.linuxgccrelease @""" + p.split("disc_out/", 1)[1] + "_args" + """

#""" + loc2 + """/source/bin/ligand_discovery.linuxgccdebug @""" + p.split("disc_out/", 1)[1] + "_args"
#for file in *; do mv "$file" `echo $file | tr ' ' '_'` ; done
#for file in *; do mv "$file" `echo $file | tr ':' '_'` ; done"
            )
            t.close()
            time.sleep(1)
            os.system("sbatch " + p.split("disc_out/", 1)[1] + ".job")
            os.chdir("../..")
        count += 1
