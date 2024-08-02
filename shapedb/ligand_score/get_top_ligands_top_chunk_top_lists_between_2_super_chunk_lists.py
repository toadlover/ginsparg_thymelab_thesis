#second step in selecting top shape-scoring ligands
#use this after the list of the top hits per super chunk are pulled by get_top_ligands_in_chunks_sub.py
#this will go through each resulting file for the superchunk (since they are generated in parallel), and look at all of the resulting files to pull the best of the best

import os,sys

#going to try to use a heap data structure to better keep the best scores
import heapq

#get args which are: the file prefix, and max number of ligands to keep


#file prefix i.e. suvo, OxB_7_shifted; do not have trailing underscore unless there is one (full file name is prefix_NN_subchunk_chunk.tar.gz)
file_prefix = str(sys.argv[1])

#value to hold the number of ligands to keep (ideally the same number as the max number in the files being looked at, but allowing the option to change if desired)
#currently using 8,000,000 when looking at 8 truncated shapes of OxA/OxB
max_ligands_to_keep = int(sys.argv[2])

#chunk addresses to look at
#example if raw 44300_44399 and 44400_44499
#later steps may have different values (not directly referencing chunk addresses)
superchunk_address_1 = str(sys.argv[3])
superchunk_address_2 = str(sys.argv[4])

#new identifiable address to write the new file as
write_address = str(sys.argv[5])

#create the heap
#each entry will be a tuple of 4 entries: shapedb score, ligand name (with conformer number), chunk, subchunk
conformer_list = []

#value to hold the currently smallest value of the conformer list so that we don't have to keep calling it every time we want it if the value has not changed
#default to -2, because the value should never get that low
smallest_value = -2

#write_file = open(str(file_prefix) + "_best_" + str(max_ligands_to_keep) + "_chunks_" + str(min_chunk_str) + "_" + str(max_chunk_str) + ".txt","w")

file1 = str(file_prefix) + "_best_" + str(max_ligands_to_keep) + "_chunks_" + str(superchunk_address_1) + ".txt"
file2 = str(file_prefix) + "_best_" + str(max_ligands_to_keep) + "_chunks_" + str(superchunk_address_2) + ".txt"
mergefile = str(file_prefix) + "_best_" + str(max_ligands_to_keep) + "_chunks_" + str(write_address) + ".txt"

#read file 1 and add to conformer list
#look at this file
print(file1)
read_file = open(file1, "r")

for line in read_file.readlines():
	#break up the line and add it to the conformer_list
	conf_score = float(line.split(",")[0])
	conf_name = line.split(",")[1]
	chunk_str = line.split(",")[2]
	subchunk_str = line.split(",")[3].strip()

	conformer_list.append([conf_score, conf_name , chunk_str, subchunk_str])

print(len(conformer_list))
read_file.close()
#read file 2 and add to conformer list
print(file2)
read_file = open(file2, "r")

for line in read_file.readlines():
	#break up the line and add it to the conformer_list
	conf_score = float(line.split(",")[0])
	conf_name = line.split(",")[1]
	chunk_str = line.split(",")[2]
	subchunk_str = line.split(",")[3].strip()

	conformer_list.append([conf_score, conf_name , chunk_str, subchunk_str])
print(len(conformer_list))

read_file.close()
#heapify and keep nlargest
conformer_heap = conformer_list

print("Heapifying list after end of current file")
heapq.heapify(conformer_heap)

#use the nlargest command to get max number back as a list to use for the next file
conformer_list = heapq.nlargest(max_ligands_to_keep,conformer_heap)
#write to new file
#once done looking at all files, write the list to a final file that contains the top ligands for this shape
write_file = open(mergefile,"w")


for i in conformer_list:
	write_file.write(str(i[0]) + "," + str(i[1]) + "," + str(i[2]) + "," + str(i[3]) + "\n")

write_file.close()
