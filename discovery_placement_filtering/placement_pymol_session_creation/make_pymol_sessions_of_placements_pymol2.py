#note, this script was initially written by ChatGPT to write most of the pymol logic, and further modified by me (Ari) to enhance use of inputs

#usage:
#python make_pymol_sessions_of_placements_pymol2.py /path/to/placements/directory/ ligandresidueindex list+string+of+residue+indices
#python make_pymol_sessions_of_placements_pymol2.py /scratch/abgvg9/discovery_results/top_1000_placement/agonist_12M_passing_placements/0 282 227+86+253+63+257

#this uses pymol2, which may be better
import os,sys
import pymol2

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
ligand_selection = "resi " + str(ligand_residue)  # Adjust as necessary for your files

# Initialize PyMOL in headless mode
with pymol2.PyMOL() as pymol:
    # Loop through each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):  # Ensure only PDB files are processed
            filepath = os.path.join(directory, filename)
            
            # Load the file
            pymol.cmd.load(filepath)
            
            # Get the object name (assumed to be the filename without the extension)
            object_name = os.path.splitext(filename)[0]
            print(object_name)
            
            # Color the entire protein
            pymol.cmd.color(protein_color, object_name)

            
            # Find and show residues within a certain distance of the ligand
            neighboring_residues = f'(byres {ligand_selection} around {distance_threshold}) and {object_name}'
            pymol.cmd.show('sticks', neighboring_residues)
            #pymol.cmd.color('elem', neighboring_residues)
            pymol.cmd.color('yellow', f'{neighboring_residues} and elem C')

            # Show the specific residues as sticks and color them
            pymol.cmd.show('sticks', f'{object_name} and resi {highlight_residues}')
            #pymol.cmd.color('elem', f'{object_name} and resi {highlight_residues}')
            pymol.cmd.color('orange', f'{object_name} and resi {highlight_residues} and elem C')

            #color the ligands magenta
            pymol.cmd.color('white', f'{ligand_selection} and elem C')

            #fixing element coloring:
            pymol.cmd.color('white', 'elem H')
            pymol.cmd.color('red', 'elem O')
            pymol.cmd.color('blue', 'elem N')
            pymol.cmd.color('green', 'elem Cl')
            pymol.cmd.color('cyan', 'elem F')
            pymol.cmd.color('brown', 'elem Br')
            pymol.cmd.color('purple', 'elem I')
            pymol.cmd.color('yellow', 'elem S')

            # Display hydrogen bonds
            pymol.cmd.dist('hbonds_temp', ligand_selection, neighboring_residues, cutoff=3.5, mode=2)
            pymol.cmd.set('dash_color', 'green', 'hbonds_temp')
            pymol.cmd.set('dash_width', 2.0)
            # Merge the hydrogen bonds into the protein-ligand object
            try:
                pymol.cmd.create(f'{object_name}', f'{object_name} or hbonds_temp')
                pymol.cmd.delete('hbonds_temp')

    # Save the session for all proteins
    pymol.cmd.save('all_proteins_session.pse')  # Save the PyMOL session