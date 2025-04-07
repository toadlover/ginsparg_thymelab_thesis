#the purpose of this script is to take in two list files of smiles and determine the overlap of molecules present in each file
#for each file, there should only be 1 smiles string per line
#this script outputs a file 1 ligand smiles per line that does overlap between the 2 files, using the smiles string as it appears in the first file
#files that are unique to a single file are not listed (this could be done in another script if we want it)
#to get this small script written quickly, logic was written with chatgpt

import os,sys
from rdkit import Chem

file1 = sys.argv[1]

file2 = sys.argv[2]


def load_and_clean_smiles(file_path):
    cleaned_smiles = set()
    with open(file_path, 'r') as f:
        for line in f:
            smi = line.strip()
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue  # skip invalid lines
            Chem.RemoveStereochemistry(mol)
            cleaned = Chem.MolToSmiles(mol, isomericSmiles=False)
            cleaned_smiles.add(cleaned)
    return cleaned_smiles

def find_common_smiles(file1, file2, output_file='matched_smiles.smi'):
    smiles1 = load_and_clean_smiles(file1)
    smiles2 = load_and_clean_smiles(file2)
    common = smiles1.intersection(smiles2)

    with open(output_file, 'w') as f:
        for smi in sorted(common):
            f.write(smi + '\n')
    print(f"Found {len(common)} matching SMILES (stereochemistry removed). Saved to {output_file}")

# Example usage
find_common_smiles(file1, file2)