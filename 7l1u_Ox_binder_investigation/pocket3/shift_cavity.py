from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def get_center(mol):
    conf = mol.GetConformer()
    coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    arr = np.array([(p.x, p.y, p.z) for p in coords])
    return arr.mean(axis=0)

# Load molecules
drug = Chem.MolFromMol2File("../../shapedb/shapedb_inputs/antagonists/empa_shifted.mol2", removeHs=False)
cavity = Chem.MolFromPDBFile("7l1u_receptor_and_orexin_empty_pocket_3_empty_space.pdb", removeHs=False, sanitize=False)

# Check load success
if drug is None:
    raise ValueError("Failed to load drug.mol2")
if cavity is None:
    raise ValueError("Failed to load cavity.pdb")

# Calculate centers
center_drug = get_center(drug)
center_cavity = get_center(cavity)

# Compute translation vector
shift = center_drug - center_cavity

# Apply translation to cavity
conf = cavity.GetConformer()
for i in range(cavity.GetNumAtoms()):
    pos = conf.GetAtomPosition(i)
    new_pos = pos.x + shift[0], pos.y + shift[1], pos.z + shift[2]
    conf.SetAtomPosition(i, new_pos)

# Save new aligned cavity
Chem.MolToPDBFile(cavity, "7l1u_receptor_and_orexin_empty_pocket_3_empty_space_aligned.pdb")

