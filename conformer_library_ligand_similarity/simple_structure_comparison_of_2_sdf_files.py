#this script serves as a simple test to get the similarity of 2 ligands using rdkit to get fingerprints and tanimoto score
#this file takes in the base ligand(s) as a sdf (only ligand in the sdf) and extracts the smiles string
#the file then takes in a second sdf file and outputs the ligand name, smiles string, and similarity score to the target

#imports
import os,sys
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

#get reference file
reference_ligands_file = sys.argv[1]

#get the other ligands list (such as split_new_named_0.sdf)
compare_ligands_list_file = sys.argv[2]

#add optional argument for using chirality in tanimoto score determination
chirality_usage = False
if len(sys.argv) >= 4:
	if sys.argv[3] == "useChirality":
		chirality_usage = True

#make a dictionary to hold the reference ligand(s) with the key as the ligand name (pulled from sdf) and value as the smiles string
reference_ligands = {}

#make a disctionary to hold the compare ligand(s) with the key as the ligand name (pulled from sdf) and value as the smiles string
compare_ligands = {}

#derive all ligand smiles and put them in the corresponding dictionaries
#reference

#get the supplier of the sdf file
supplier = Chem.SDMolSupplier(reference_ligands_file)
#iterate over the molecules in the supplier
for mol in supplier:
    #skip invalid molecules
    if mol is not None:
        #derive the smiles string
        smiles = Chem.MolToSmiles(mol)
        # Extract the name or identifier (usually stored in the "MOL" block or as a property)
        lig_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unnamed_Ligand"
        
        #add the reference ligand to the dictionary
        reference_ligands[lig_name] = smiles
print("Got reference ligands.")

#repeat for compare ligands file
#get the supplier of the sdf file
supplier = Chem.SDMolSupplier(compare_ligands_list_file)
#iterate over the molecules in the supplier
for mol in supplier:
    #skip invalid molecules
    if mol is not None:
        #derive the smiles string
        smiles = Chem.MolToSmiles(mol)
        # Extract the name or identifier (usually stored in the "MOL" block or as a property)
        lig_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Unnamed_Ligand"
        
        #add the reference ligand to the dictionary
        compare_ligands[lig_name] = smiles
print("Got compare ligands.")

#iterate over each reference ligand
#for each ligand, prepare a csv file that is named after the reference ligand, and then lists each compare ligand
for curr_ref_lig in reference_ligands.keys():
	
	print("On reference ligand " + curr_ref_lig + ": " + reference_ligands[curr_ref_lig])

	#make the output csv
	write_file = open(curr_ref_lig + "_smiles_similarity.csv","w")

	#push smiles back into a mol object
	ref_mol_smiles = Chem.MolFromSmiles(reference_ligands[curr_ref_lig])

	#get the fingerprints of the reference
	ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol_smiles, radius=2, useChirality = chirality_usage)

	#iterate over each compare ligand
	for curr_comp_lig in compare_ligands.keys():
		#compare the current compar ligand against the current reference ligand

		#push smiles back into a mol object
		comp_mol_smiles = Chem.MolFromSmiles(compare_ligands[curr_comp_lig])

		#get the compare ligand fingerprints
		comp_fp = AllChem.GetMorganFingerprintAsBitVect(comp_mol_smiles, radius=2, useChirality = chirality_usage)
	    
		#get the tanimoto similarity by fingerprints
		similarity = DataStructs.TanimotoSimilarity(ref_fp, comp_fp)

		#write the data to the csv
		write_file.write(curr_comp_lig + "," + str(similarity) + "," + compare_ligands[curr_comp_lig] + "\n")