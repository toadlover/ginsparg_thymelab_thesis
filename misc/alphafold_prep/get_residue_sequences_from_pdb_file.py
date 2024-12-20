#this script takes in a single pdb file, and prints a residue sequence of each chain.
#note, most of this code was written by ChatGPT as a quick way to get this logic, and was then adapted to suit my needs
from rdkit import Chem
from rdkit.Chem import rdmolfiles

def extract_chains_from_pdb(pdb_file):
    mol = Chem.MolFromPDBFile(pdb_file)
    if mol is None:
        print("Failed to load the PDB file.")
        return None
    
    chains = {}
    for atom in mol.GetAtoms():
        # Get chain ID from atom properties
        chain_id = atom.GetPDBResidueInfo().GetChainId()
        if chain_id not in chains:
            chains[chain_id] = Chem.RWMol()
        chains[chain_id].AddAtom(atom)
    
    # Convert to read-only molecules
    chain_mols = {chain_id: chain_mol.GetMol() for chain_id, chain_mol in chains.items()}
    return chain_mols

# Example usage
pdb_file = sys.argv[1]
chains = extract_chains_from_pdb(pdb_file)

if chains:
    for chain_id, chain_mol in chains.items():
        print(f"Chain {chain_id}:")
        print(Chem.MolToSmiles(chain_mol))
