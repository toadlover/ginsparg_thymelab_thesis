import os, sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

smiles_strings = {}

for i in range(len(sys.argv)):
    if i >= 1:
        smiles_strings[i] = sys.argv[i]

# Function to get molecular formula and mass
def smiles_to_formula_and_mass(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    if mol:
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mass = rdMolDescriptors.CalcExactMolWt(mol)
        return formula, mass
    return None, None

# Calculate formulas and masses
results = {}
for name, smiles in smiles_strings.items():
    formula, mass = smiles_to_formula_and_mass(smiles)
    results[name] = {"SMILES": smiles, "Formula": formula, "Exact Mass": mass}

# Pretty print the results
for idx, data in results.items():
    print(f"Input {idx}:")
    print(f"  SMILES: {data['SMILES']}")
    print(f"  Formula: {data['Formula']}")
    print(f"  Exact Mass: {data['Exact Mass']}\n")