import os

#Need to add path to unpacked "all" folder
p = "conversion_test_space"

dirlist = [ item for item in os.listdir(p) if os.path.isdir(os.path.join(p, item)) ]

os.chdir(p)
for d in dirlist:
    os.chdir(d)
    os.system("mv receptor.pdb " + d + ".pdb")
    n = open(d + "_args", "x")
    n.write(
"""#optional flag to output generated params file as a pdb (to confirm that the params file was converted correctly)
-write_param_to_pdb true
#optional flag to fully write a custom name of the ligand into the params file (helps make sure fill name isn't lost in  pipeline)
-ligand_full_name ''' + d + "_lig" + '''
#optional flag to custom override the 1 and 3 letter codes for the ligand in the params file
#-ligand_1_letter_code N
#-ligand_3_letter_code New
#indicates if ligand  is an  amino  acid  or not (will almost never be used)
#-mol2_amino_acid ALA
#MOST IMPORTANT FLAG! Source mol2 file with path
-mol2_file crystal_ligand.mol2
-params_custom_file_name """ + d + "_lig.params"
    )
    
    job_name = d + ".job"
    job = open(job_name, "x")
    job.write("#!/bin/bash\n")
    job.write("#SBATCH -p express # Partition to submit to\n")
    job.write("#SBATCH -n 1 # Number of cores requested\n")
    job.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
    job.write("#SBATCH -t 120 # Runtime in minutes\n")
    job.write("#SBATCH --mem=20000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
    job.write("#SBATCH -o " + d +".out # Standard out goes to this file\n")
    job.write("#SBATCH -e " + d +".err # Standard err goes to this filehostname\n")
    job.write("/data/project/thymelab/running_Rosetta/ari_work/Rosetta_Code_copy/main/source/bin/"
              "mol2_to_params.linuxgccrelease @" + d + "_args" + "\n")
    job.close()
    os.system("sbatch " + job_name)
    os.chdir("..")



