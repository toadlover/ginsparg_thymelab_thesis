#note, this script was initially written by ChatGPT to write most of the pymol logic, and further modified by me (Ari) to enhance use of inputs

#usage:
#python make_pymol_sessions_of_placements.py /path/to/placements/directory/ ligandresidueindex list+string+of+residue+indices
#python make_pymol_sessions_of_placements.py /scratch/abgvg9/discovery_results/top_1000_placement/agonist_12M_passing_placements/0 282 227+86+253+63+257

import os,sys
import pymol
from pymol import cmd

# Initialize PyMOL
pymol.finish_launching()

# Define the directory containing your files as a command line argument
#directory = 'path/to/your/files'
directory = sys.argv[1]

#have the directory end with a backslash if it doesn't from the input
if directory.endswith("/") == False:
    directory = directory + "/"

#ligand residue index of interest
ligand_residue = sys.argv[2]

# Define the specific residues you want to highlight (e.g., 86, 227, 342)
highlight_residues = sys.argv[3]

# Define the color for the entire protein
protein_color = 'cyan'

# Define the distance threshold (in angstroms)
distance_threshold = 5.0

# Define the selection around which to find neighboring residues (e.g., ligand)
ligand_selection = 'ligand'  # Adjust as necessary for your files

# Loop through each file in the directory
for filename in os.listdir(directory):
    if filename.endswith('.pdb'):  # Ensure only PDB files are processed
        filepath = os.path.join(directory, filename)
        
        # Load the file
        cmd.load(filepath)
        
        # Get the object name (assumed to be the filename without the extension)
        object_name = os.path.splitext(filename)[0]
        
        # Color the entire protein
        cmd.color(protein_color, object_name)
        
        # Find and show residues within a certain distance of the ligand
        neighboring_residues = f'(byres {ligand_selection} around {distance_threshold}) and {object_name}'
        cmd.show('sticks', neighboring_residues)
        cmd.color('yellow', neighboring_residues)

        # Show the specific residues as sticks and color them
        cmd.show('sticks', f'{object_name} and resi {highlight_residues}')
        cmd.color('red', f'{object_name} and resi {highlight_residues}')
        


# Save the session for all proteins
cmd.save(directory + 'all_proteins_session.pse')  # Save the PyMOL session

# Finish PyMOL
pymol.cmd.quit()