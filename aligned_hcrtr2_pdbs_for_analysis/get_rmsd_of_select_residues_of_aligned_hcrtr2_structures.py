#the purpose of this script is to take the coordinates of select binding pocket residues and calculate the rmsd of the residue against 4s0v (after being previously aligned to 4s0v)

#this is intended to be a 1-off script to process data, but could be adapted in the future
#this script will produce a csv file that I can make a heatmap off of to determine the rmsd per residue per structure against 4s0v

#import packages
import os,sys

#read in the list of residues of interest from 4s0v
#store to a list of tuples where index 0 is the residue 3 letter code and index 1 is the residue index
residues_of_interest = []

#read in the list
pocket_residue_file = open("4s0v_pocket_residues_unique.txt", "r")
for line in pocket_residue_file.readlines():
	#strip the newline
	stripped_line = line.strip()
	#make sure we are looking at a good line
	if len(stripped_line.split()) == 4:
		#then assign

		residue = stripped_line.split()[2]
		index = stripped_line.split()[3]
		residues_of_interest = [residue,index]

#sanity check print of the residues of interest
for index in residues_of_interest:
	print(index)