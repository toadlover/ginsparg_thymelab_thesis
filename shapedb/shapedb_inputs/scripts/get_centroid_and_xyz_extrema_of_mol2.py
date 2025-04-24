#This script serves to determine the centroid coordinate and xyz extrema of an inputted molecule (defines a rectangular prism around it)
#the rectangular coordinate extrema are defined as 8 coordinates (probably less useful) and as side lengths for the prism
#currently only works for mol2 files (with a single molecule), but could probably work for any file type if I worked babel into it; could slightly modify the script for sdf and pdb versions if needed

import os,sys

#get the file name
mol2_file = sys.argv[1]

#bool control variable to determine when we are on atom lines so we take in coordinates
reading_atoms = False

#variables to hold xyz extreme and centroid
#define centroid at end after we get extrema
controid = []
#use 'x' as a placeholder to trigger initial seeding
xmin = 'x'
xmax = 'x'
ymin = 'x'
ymax = 'x'
zmin = 'x'
zmax = 'x'

dummy_str = 'x'

#read the file
read_file = open(mol2_file,"r")

for line in read_file.readlines():
	#determine when we are on atom lines (with coordinates)
	
	#if we hit any other @<TRIPOS>line, it is not atoms (unless it is the ATOM line, in which we flip to true)
	if "@<TRIPOS>" in line:
		reading_atoms = False	

	#read the atom lines
	if reading_atoms:
		#parse the line to get the data
		x = float(line.split()[2])
		y = float(line.split()[3])
		z = float(line.split()[4])

		#determine if any coordinates correspond to an extrema
		
		#automatically assign to current value if it is the first value encountered
		if type(dummy_str) == type(xmin):
			xmin = x
		if type(dummy_str) == type(xmax):
			xmax = x
		if type(dummy_str) == type(ymin):
			ymin = y
		if type(dummy_str) == type(ymax):
			ymax = y
		if type(dummy_str) == type(zmin):
			zmin = z
		if type(dummy_str) == type(zmax):
			zmax = z

		if x < xmin:
			xmin = x

		if x > xmax:
			xmax = x

		if y < ymin:
			ymin = y

		if y > ymax:
			ymax = y

		if z < zmin:
			zmin = z

		if z > zmax:
			zmax = z

	#if we hit the @<TRIPOS>ATOM line, we can start reading lines after as atom lines
	#keep control at end since the atoms start in the next line
	if "@<TRIPOS>ATOM" in line:
		reading_atoms = True

#determine the side lengths
xside = xmax - xmin
yside = ymax - ymin
zside = zmax - zmin

#get half ot the side lengths and convert to int
xrad = int(xside/2)
yrad = int(yside/2)
zrad = int(zside/2)

#determine the centroid, cast to int
#centroid coordinate defined as max - 0.5 * side length
centroid = [int(xmax-(xside/2)),int(ymax-(yside/2)),int(zmax-(zside/2))]

#print out the centroid and sidelengths in the args format to easily add to a rosetta args file
print("-binding_pocket_center_sf " + str(centroid[0]) + "," + str(centroid[1]) + "," + str(centroid[2]))
print("-binding_pocket_dimensions_sf " + str(xrad) + "," + str(yrad) + "," + str(zrad))

