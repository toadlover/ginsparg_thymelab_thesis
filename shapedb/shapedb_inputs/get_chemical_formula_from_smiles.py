#note, this script was written with the help of chatgpt as a simple tool to take in a SMILES string and get the chemical formula of the string

import os,sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

smiles_strings = {}

for i in range(len(sys.argv)):
    if i >= 1:
        smiles_strings[i] = sys.argv[i]

# Function to get molecular formula
def smiles_to_formula(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    if mol:
        return rdMolDescriptors.CalcMolFormula(mol)
    return None

# Calculate formulas
formulas = {name: smiles_to_formula(smiles) for name, smiles in smiles_strings.items()}
print(formulas)

