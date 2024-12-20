import os,sys

#use Chem in rdkit to derive SMILES string of placed ligands
from rdkit import Chem

#get ligand file (and path) that you want to derive smiles for
#smiles file that is written will take the name up to the first period
#for the sake of simplicity, this only works for mol2 and sdf ligands
input_file = sys.argv[1]

#sanity check for file type
if input_file.endswith(".sdf") == False and input_file.endswith(".mol2") == False and input_file.endswith(".pdb") == False:
	print("File input of " + input_file + " is not of type pdb, sdf, or mol2, and we will not mess with it")
	quit()

file_prefix = input_file.split("/")[len(input_file.split("/")) - 1].split(".")[0]

#optional, location to write the smiles file to. otherwise default to writing to current location
output_path = ""
if len(sys.argv) >= 3:
	output_path = sys.argv[2]

	#append a slash to output_path if it does not end with one
	if output_path.endswith("/") == False:
		output_path = output_path + "/"

#optional, replace the file prefix with custom as 3rd argument
if len(sys.argv) >= 4:
	file_prefix = sys.argv[3]

#derive the smiles string

mol = ""

if input_file.endswith(".sdf"):
	mol = Chem.SDMolSupplier(input_file)
if input_file.endswith(".mol2"):
	mol = Chem.MolFromMol2File(input_file)
if input_file.endswith(".pdb"):
	mol = Chem.MolFromPDBFile(input_file)


if mol is not None:  # Check for valid molecule
	smiles = Chem.MolToSmiles(mol)

	#write the smiles to a new file
	#print(output_path + file_prefix + ".smi",smiles)
	write_file = open(output_path + file_prefix + ".smi","w")
	write_file.write(smiles)
	write_file.close()