#this is a small helper script that reads in an input text file that has 1 file per line (ideally with a leading path before the file), and copies it to a given target location

#imports
import os,sys

in_file = sys.argv[1]
location = sys.argv[2]

#read the file
read_file = open(in_file,"r")

for line in read_file.readlines():
	#strip the line
	stripped_line = line.strip()

	#copy the file to the target location
	os.system("cp " + stripped_line + " " + location)