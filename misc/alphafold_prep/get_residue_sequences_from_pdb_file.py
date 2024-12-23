#this script takes in a single pdb file, and prints a residue sequence of each chain.
#note, most of this code was written by ChatGPT as a quick way to get this logic, and was then adapted to suit my needs

import os,sys

#this is a basic script that just reads the raw file, looks at ATOM lines, and determines the residue sequence based on the order that the residues appear
#each chain is individually listed, as alphafold wants that

pdb_file = sys.argv[1]

#dictionary composed by chatgpt to help handle conversion of all cacnonical and some noncanonical residues from 3 letter code to 1 letter code
amino_acid_codes = {
    # Canonical amino acids
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    
    # Common noncanonical amino acids
    "ASX": "B",  # Aspartic acid or asparagine
    "GLX": "Z",  # Glutamic acid or glutamine
    "SEC": "U",  # Selenocysteine
    "PYL": "O",  # Pyrrolysine
    "UNK": "X",  # Unknown amino acid

    # Other nonstandard or modified residues
    "MSE": "M",  # Selenomethionine
    "CME": "C",  # Carboxymethyl cysteine
    "CSS": "C",  # Thioether cysteine
    "CSO": "C",  # Cysteine-sulfinic acid
    "HIP": "H",  # Protonated histidine
    "HIE": "H",  # Neutral histidine (epsilon-tautomer)
    "HID": "H",  # Neutral histidine (delta-tautomer)
    "TPO": "T",  # Phosphothreonine
    "SEP": "S",  # Phosphoserine
    "PTR": "Y",  # Phosphotyrosine
    "LEU": "L",  # Leucine
    "CYX": "C",  # Disulfide-linked cysteine
    "PCA": "E",  # Pyroglutamate
}

#make a dictionary of chains, the key is the chain id (from some example pdb files I have, a space " " can be a chain id)
#the values for each key are a tuple that contains the residue indices (for ensuring duplicates are not recorded) and the 1 letter code sequence of residues
chains = {}

#read through the pdb file and get the data
read_file = open(pdb_file, "r")

for line in read_file.readlines():
    #only look at lines starting with "ATOM"
    if line.startswith("ATOM"):
        print(line)

        #extract the chain, 3 letter code, and residue index
        chain_id = line[21]
        residue_3_code = line[17] + line[18] + line[19]
        residue_index = int(str(line[23] + line[24] + line[25]).strip())

        print(chain_id,residue_3_code,residue_index)

        #determine how/if to add the residue 1 letter code to the chain data

        #initial check to see if we have data on this chain (or start a new chain)

        #check if we have already seen this index

        #

