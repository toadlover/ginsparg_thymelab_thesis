import os,sys

file = sys.argv[1]
compare_sdf = sys.argv[2]
#extract the file name for the output scores
filebase = file.split("/")[len(file.split("/")) - 1]

dir_and_sub = ""

if len(file.split("/")[len(file.split("/")) - 1].split("_")) == 2:
	dir_and_sub = str(filebase.split("_")[0]) + "_" + str(filebase.split("_")[1])
else:
	dir_and_sub = str(filebase.split(".")[0])

#extract the base filename from the compare sdf
sdfbase = compare_sdf.split("/")[len(compare_sdf.split("/")) - 1]
sdf_base = sdfbase.split(".")[0]

#read in the conformer list text file
readfile = open(file, "r")

line_counter = 0
for line in readfile.readlines():
	line_counter = line_counter + 1

#print(file, compare_sdf, dir_and_sub, sdf_base)

#run nnsearch
#db.db should be located where you are
os.system("/pharmit/src/build/shapedb -NNSearch -k " + str(line_counter) + " -ligand " + compare_sdf + " -db db.db -print > " + dir_and_sub + "_scored_confs_against_" + sdf_base + ".txt")