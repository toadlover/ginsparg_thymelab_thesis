import os,sys

#use Chem in rdkit to derive SMILES string of placed ligands
from rdkit import Chem

#get the directory to runconformator on
working_location = sys.argv[1]

#get the location of the conformator executable
#include path and executable (in casr the executable name is different for some reason)
conformator_executable = sys.argv[2]

#move to working location
os.chdir(working_location)

#run through each sdf file in the working location
for r,d,f in os.walk(working_location):
	for file in f:
		if file.endswith(".sdf") and r == working_location:
			
			#tracer
			print(file)

			#break up the file name to remove extension
			file_prefix = file.split(".")[0]

			#derive the msiles string
			supplier = Chem.SDMolSupplier(sdf_file)
			for mol in supplier:
				if mol is not None:  # Check for valid molecule
					smiles = Chem.MolToSmiles(mol)

					#write the smiles to a new file
					write_file = open(file_prefix + ".smi","w")
					write_file.write(smiles)
					write_file.close()

			#run conformator on the smiles file
			os.system(conformator_executable + " -i " + file_prefix + ".smi" + " -o " + file_prefix + "_confs.sdf --keep3d --hydrogens -v 0")
#			os.system("/conformator_for_container/conformator_1.2.1/conformator -i " + starting_file + " -o confs.sdf --keep3d --hydrogens -n 15 -v 0")