#the purpose of this script is to be a bypass of some slower methods
#it is intended to be used after running ginsparg_thymelab_thesis/ discovery_placement_filtering/placement_pymol_session_creation/convert_placement_pdb_ligands_to_sdf.py and be used on a directory of sdf files
#this script will make a dummy csv file containing a list of smiles with ligand names to then be fed into ginsparg_thymelab_thesis/conformer_library_ligand_similarity/prepare_conformers_and_params_from_smiles_comparison_csv.py and reduce steps (useful for pushing the data up to a bucket and run on osg)

#package imports
from rdkit import Chem
import os,sys

#get the folder of sdf files as command line argument
sdf_folder = sys.argv[1]

#if folder doesn't end with a /, make it end with one
if sdf_folder.endswith("/") == False:
	sdf_folder = sdf_folder + "/"


#declare list of smiles, will be a tuple of the ligand name and the smiles
smiles_list = []

#iterate over the sdf files in the folder to get the smiles string
for file in os.listdir(sdf_folder):
    if file.endswith(".sdf"):
        sdf_path = os.path.join(sdf_folder, file)
        suppl = Chem.SDMolSupplier(sdf_path)
        
        for mol in suppl:
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                
                #derive the ligand name as whatever is before the first period and underscore
                #do not pass on underscores for 64b library ligands
                lig_name = file.split(".")[0]

                smiles_list.append([lig_name,smiles])
                break  # Only extract the first (single) ligand

#now write a dummy list to the location
write_file = open(sdf_folder + "dummy_list.csv", "w")

#write dummy csv lines that have the ligand name and smiles string
for lig in smiles_list:
	#write ligand name to 3rd entry and smiles to 5th
	write_file.write(",," + lig[0] + ",," + lig[1] + "\n")

