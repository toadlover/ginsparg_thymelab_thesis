#works on all pdb files in the current directory, and appends "lig_" to the front of the source name
for f in *.pdb; do grep '^HETATM' "$f" | grep -Ev ' HOH| WAT' > "lig_${f}"; done
