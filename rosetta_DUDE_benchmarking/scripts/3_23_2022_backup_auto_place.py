import os
import re
import os.path
import sys

loc1 = str(sys.argv[1])
loc2 = str(sys.argv[2])
fa_rep_cutoff = str(sys.argv[3])
fa_atr_cutoff = str(sys.argv[4])

os.mkdir("disc_out")
for f in os.listdir("/data/user/abgvg9/benchmark/motif_collection/files/neighbor_out"):
    d = f.replace("output_", "")
    d = d.replace(".txt", "")
    o = open("/data/user/abgvg9/benchmark/motif_collection/files/neighbor_out/" + f, "r")
    count = 0
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

-motif_filename """ + loc1 + """

#CHANGE THIS
-protein_discovery_locus """ + n + """

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
loc2 + """/main/source/bin/ligand_discovery.linuxgccrelease @""" + p.split("disc_out/", 1)[1] + "_args" + """
#""" + loc2 + """/main/source/bin/ligand_discovery.linuxgccdebug @""" + p.split("disc_out/", 1)[1] + "_args"
            )
            t.close()
            os.system("sbatch " + p.split("disc_out/", 1)[1] + ".job")
            os.chdir("../..")
        count += 1
