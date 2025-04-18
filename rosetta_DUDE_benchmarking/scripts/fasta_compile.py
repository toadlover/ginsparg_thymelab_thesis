import os

filenames = []
for f in os.listdir("fasta"):
    filenames.append(f)
with open('fasta/fasta_compiled.fasta', 'w') as outfile:
    for fname in filenames:
        with open("fasta/" + fname) as infile:
            outfile.write(infile.read())
