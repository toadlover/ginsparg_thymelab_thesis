#the purpose of this script is to run the smiles comparison script on a superchunk of the library
#the script will take in the superchunk to process (0-530); 0 corresponds to chunks 00000-00099, 1 corresponds to 00100-00199, ... , 530 corresponds to 53000-53084 (53084 is the highest chunk)
#the script will also take in arguments for the comparison sdf file
#the script also takes in a number for the top X most similar ligands to isolate from this data

#import packages
import os,sys
#import heapq to heapify the most similar ligands, which we can track and output for later
import heapq

#get the superchunk
superchunk = int(sys.argv[1])

#get the reference ligand sdf
#give the absolute path to the ligand file!
reference_ligand = sys.argv[2]

#derive the reference ligand name for later use
reference_ligand_name = ""
reference_ligand_file = open(reference_ligand,"r")
for line in reference_ligand_file.readlines():
	#the ligand name should be on line 1, so break after
	reference_ligand_name = line.strip()
	break
reference_ligand_file.close()

#get the number for the top X ligands to consider keeping
top_ligand_amount = int(sys.argv[3])

#get the path of this script so that the comparison script can be called
script_dir = os.path.dirname(os.path.abspath(__file__))

#declare the compare script is a string for when it gets called
compare_script = script_dir + "/simple_structure_comparison_of_2_sdf_files.py"

#declare an empty list to hold the top x ligands
top_ligands = []

#initial behavior, create a directory for the superchunk to help with organization
os.system("mkdir " + str(superchunk))

#move into the directory
os.chdir(str(superchunk))

#derive the minimum and maximum chunk based on the superchunk
#min is superchunk * 100, max is min + 100 (will drop by 1 when using the range operator)
min_chunk = superchunk * 100
#max_chunk = min_chunk + 100
max_chunk = min_chunk + 1

#now, iterate from 0-99 and iterate over each chunk and get the similarity of each ligand in each subchunk from each chunk
for i in range(min_chunk,max_chunk):
	#create a directory for the corresponding chunk
	#identify the chunk ID and write to a string
	chunk = str(i)

	#add leading zeroes until chunk string is 5 characters long
	while len(chunk) < 5:
		chunk = "0" + chunk

	print("On chunk " + chunk)

	#create a directory for the chunk
	os.system("mkdir " + chunk)

	#move into the chunk
	os.chdir(chunk)

	#iterate over the chunk for each subchunk from 0-9
	for j in range(0,10):
		#attempt to pull the subchunk down
		os.system("s3cmd get s3://ariosg/ligand_library/" + chunk + "/for_s3/split_new_named_" + str(j) + ".sdf.tar.gz")

		#unzip the file
		os.system("tar -xzf split_new_named_" + str(j) + ".sdf.tar.gz")

		#run the comparison script on the reference file and the newly acquired file
		os.system("python " + compare_script + " " + reference_ligand + " split_new_named_" + str(j) + ".sdf")

		#we now have the comparison csv data, read it and rewrite the file
		compare_file_in = open(reference_ligand_name + "_smiles_similarity.csv", "r")
		#also open a write stream to rename the file and write the data with the chunk and subchunk
		compare_file_out = open(str(j) + "_" + reference_ligand_name + "_smiles_similarity.csv", "w")

		for line in compare_file_in.readlines():
			
			#break if the line does not have enough commas
			if len(line.split(",")) < 3:
				break

			#break up the line into a tuple
			#order of elements will now be the chunk,subchunk,ligname,similarity_score,smiles_string
			#declare temporary tuple to work with
			line_tuple = []

			#write the chunk and subchunk to the tuple
			line_tuple.append(str(chunk))
			line_tuple.append(str(j))

			#create stripped line
			stripped_line = line.strip()

			#append the ligand name, similarity score (take the negative for use with heap), and smiles string
			line_tuple.append(str(stripped_line.split(",")[0]))
			line_tuple.append(str(float(stripped_line.split(",")[1])))
			line_tuple.append(str(stripped_line.split(",")[2]))

			#we now have the line tuple, attempt to integrate it into the top x heap
			heapq.heappush(top_ligands, (float(stripped_line.split(",")[1]) * -1 ,line_tuple))

			#if the length of the heap is above the top x, pop off the worst element
			if len(top_ligands) > top_ligand_amount:
				heapq.heappop(top_ligands)

			#write the revised line to the new csv
			compare_file_out.write(str(chunk) + "," + str(j) + "," + str(stripped_line.split(",")[0]) + "," + str(stripped_line.split(",")[1]) + "," + str(stripped_line.split(",")[2]) + "\n")

		#we are now done processing the original csv, delete it and the corresponding sdf files
		os.system("rm " + reference_ligand_name + "_smiles_similarity.csv *.sdf*")

	#we are now done with all subchunks
	#temporary test print of the top 10 in the heap (will likely remove later)
	print("Current top 10 closest ligands from this search: ")
	counter = 0
	for lig in top_ligands:
		print(lig)
		counter = counter + 1
		if counter >= 10:
			break

	#move up one level so that we can move to the next chunk (or move into program end behavior)
	os.chdir("..")

#in this directory, we now have our top x most similar ligands from this superchunk (in which other subchunks may outperform and a final coalescence is needed)
#write the best from this superchunk into a new csv file
best_ligs_file = open("superchunk_top_" + str(top_ligand_amount) + ".csv", "w")

print("Done collecting all data, writing top " + str(top_ligand_amount) + " closest matching ligands from this superchunk")

for lig in top_ligands:
	best_ligs_file.write(str(lig[1][0]) + "," + str(lig[1][1]) + "," + str(lig[1][2]) + "," + str(float(lig[1][3])) + "," + str(lig[1][4]) + "\n")
