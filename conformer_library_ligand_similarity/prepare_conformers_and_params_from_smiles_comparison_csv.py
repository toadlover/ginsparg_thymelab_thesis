#the purpose of this script is to take in a csv file that has the smiles strings of ligands of interest
#we take the smiles, derive sdf files of the smiles, and use conformator to build a diverse conformer set (use default cutoff of 150), and make rosetta params of the conformers to investigate with Rosetta
#each ligand will get its own folder, and will be uploaded to a bucket location of choice

#csv line format:
#chunk,subchunk,ligname,similarityscore,smiles

#needed inputs
#results csv file
#path to conformator executable
#path to working copy of Rosetta molfile_to_params.py
#bucket location to push compressed test_params directories

#package imports
#this utilizes rdkit, as well as having conformator, Rosetta, and s3cmd
import os,sys
from rdkit import Chem
from rdkit.Chem import AllChem

#read in mandatory inputs in order

#results csv file, can include path
csv_file = sys.argv[1]

#conformator executable
#
conformator_executable = sys.argv[2]

#