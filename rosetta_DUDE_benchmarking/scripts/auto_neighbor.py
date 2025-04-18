import os

#Need to add path to unpacked "all" folder
p = "all"

dirlist = [ item for item in os.listdir(p) if os.path.isdir(os.path.join(p, item)) ]
os.system("source activate rosetta_tools_env")
os.chdir(p)
for d in dirlist:
    os.chdir(d)
    os.system("python3 /data/project/thymelab/running_Rosetta/ari_work/Rosetta_tools/tools/"
              "protein_tools/scripts/clean_pdb.py " + d + ".pdb ignorechain")
    os.system("mv " + d + ".pdb " + d + "_orig.pdb")
    os.system("mv " + d + "_ignorechain.pdb " + d + ".pdb")
    os.system("wget http://dude.docking.org/targets/" + d + "/docking/xtal-lig.pdb")
    os.system("mv xtal-lig.pdb " + d + "-lig.pdb")
    n = open("list.txt", "x")
    n.write(d + ".pdb \n" + d + "-lig.pdb")
    n.close()
    os.system("cp -r /data/user/ishaanp@sas.upenn.edu/test/test_params .")
    os.system("cp /data/user/ishaanp@sas.upenn.edu/test/neighbor_args .")
    os.chdir("test_params")
    os.system("cp ../" + d + "_lig.params .")
    n = open("residue_types.txt", "x")
    n.write(
"""## the atom_type_set and mm-atom_type_set to be used for the subsequent paramet$
TYPE_SET_MODE full_atom
ATOM_TYPE_SET fa_standard
ELEMENT_SET default
MM_ATOM_TYPE_SET fa_standard
ORBITAL_TYPE_SET fa_standard
##
## Test params files \n""" +
d + "_lig.params"
    )
    n.close()
    os.chdir("..")
    job_name = d + ".job"
    job = open(job_name, "x")
    job.write("#!/bin/bash\n")
    job.write("#SBATCH -p express # Partition to submit to\n")
    job.write("#SBATCH -n 1 # Number of cores requested\n")
    job.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
    job.write("#SBATCH -t 10 # Runtime in minutes\n")
    job.write("#SBATCH --mem=10000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
    job.write("#SBATCH -o " + d +".out # Standard out goes to this file\n")
    job.write("#SBATCH -e " + d +".err # Standard err goes to this filehostname\n")
    job.write("/data/project/thymelab/running_Rosetta/ari_work/Rosetta_Code_copy/main/source/bin/neighbor.linuxgccrelease"
              " @neighbor_args > " + d + "_out.txt")
    job.close()
    os.system("sbatch " + job_name)
    os.chdir("..")



