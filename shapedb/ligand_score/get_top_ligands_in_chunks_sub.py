import os,sys

#going to try to use a heap data structure to better keep the best scores
import heapq

#get args which are: the file prefix, chunk finder (int between 0 and 530), and max number of ligands to keep


#file prefix i.e. suvo, OxB_7_shifted; do not have trailing underscore unless there is one (full file name is prefix_NN_subchunk_chunk.tar.gz)
file_prefix = str(sys.argv[1])

#for the chunk finder, the value will be the minimum chunk, and we will look at 1000 chunks from there
chunk_finder = int(sys.argv[2])

#value to hold the number of ligands to keep
max_ligands_to_keep = int(sys.argv[3])

#define the min and max chunks to look at
#i.e. 0 -> 0, 10 -> 1000, 530 -> 53,000
min_chunk = chunk_finder * 100

#define max chunk as the min + 99
#max_chunk = min_chunk + 99
#testing statement to only look at 10 chunks for now while testing
max_chunk = min_chunk + 99

#create the heap
#each entry will be a tuple of 4 entries: shapedb score, ligand name (with conformer number), chunk, subchunk
conformer_list = []

#heapq.heapify(conformer_list)

#value to hold the currently smallest value of the conformer list so that we don't have to keep calling it every time we want it if the value has not changed
#default to -2, because the value should never get that low
smallest_value = -2

#i.e. 0 will have the range 0-99, 1 is 100-199, 10 is 1000-1099, 530 is 53,000-53,099 (this will be specifically cut to 53084 since that is the max)
#make strings for the min and max chunk
min_chunk_str = str(min_chunk)
while len(min_chunk_str) < 5:
	min_chunk_str = "0" + min_chunk_str

max_chunk_str = str(max_chunk)
while len(max_chunk_str) < 5:
	max_chunk_str = "0" + max_chunk_str

#iterate through each chunk and its subchunks of ligands to pull the best scoring ligands for the given shape file prefix
for i in range(min_chunk,max_chunk + 1):
	
	#if the current chunk (i) is greater than 53084 (highest chunk), break the loop
	if i > 53084:
		print("We are exceeding the size of the library, breaking the loop.")
		break

	#construct a 5 digit string to access chunks (lead with zeroes)
	chunk_str = str(i)

	while len(chunk_str) < 5:
		chunk_str = "0" + chunk_str

	#iterate through each subchunk (0-9)
	for j in range(0,10):
		
		subchunk_str = str(j)

		#pull the zipped file
		# s3cmd get s3://ariosg/ligand_library/00000/for_s3/OxB_7_shifted_NN_9_00000.tar.gz
		os.system("s3cmd get s3://ariosg/ligand_library/" + chunk_str + "/for_s3/" + file_prefix + "_NN_" + subchunk_str + "_" + chunk_str + ".tar.gz")

		#unzip the file
		os.system("tar -xzvf " + file_prefix + "_NN_" + subchunk_str + "_" + chunk_str + ".tar.gz")

		#open the text file to read
		try:
			read_file = open(chunk_str + "_" + subchunk_str + "_scored_confs_against_" + file_prefix + ".txt" , "r")
		except:
			print(chunk_str + "_" + subchunk_str + "_scored_confs_against_" + file_prefix + ".txt is missing")
			continue

		#read through the file
		for line in read_file.readlines():
			
			#print(line,i,j)

			#get conformer name and score
			
			#shapedb scores are multiplied by negative 1 so that the worse scoring placements are what gets popped off (since python heaps pop the lowest value)
			conf_score = float(line.split()[1].strip()) * -1
			conf_name = str(line.split()[0])

			#add entry to the heap if the size of the heap is below the size of max_ligands_to_keep
			
			#if len(conformer_list) < max_ligands_to_keep:
				#heapq.heappush(conformer_list, [conf_score, conf_name , chunk_str, subchunk_str] )

				#derive the smallest value in the heap
				#smallest_value = heapq.nsmallest(1, conformer_list)[0][0]

				#quicker and smarter approach to initially building the heap up to the size of max_ligands_to keep
				#build a a list (no need to sort every time while we build and are agnostic to the worst value)
				#heapify the list once we cap at the max desired size for the first time

			conformer_list.append([conf_score, conf_name , chunk_str, subchunk_str])

				#if len(conformer_list) == max_ligands_to_keep:
				#	#print("hit the max, time to heapify")
				#	heapq.heapify(conformer_list)
				#	#print("we are heapified")

			#if we hit the max number of entries in the heap, check if the current score would 
			#else:
			#	#print(smallest_value)
			#	#if the conformer score is better than the worst score in the list, perform a heappushpop and get the new smallest value
			#	if conf_score > smallest_value:
			#		heapq.heappushpop(conformer_list, [conf_score, conf_name , chunk_str, subchunk_str])

			#		smallest_value = heapq.nsmallest(1, conformer_list)[0][0]




		#close the file
		read_file.close()

		#remove the tar file and text file
		os.system("rm " + file_prefix + "_NN_" + subchunk_str + "_" + chunk_str + ".tar.gz " + chunk_str + "_" + subchunk_str + "_scored_confs_against_" + file_prefix + ".txt")

	#even smarter approach, only heapify at the end of the file after adding everything from the file into the list, and then pop off the lowest if we are beyond the max (and don't even heapify if we are below the max)
	#even more smart, only do this step when completing a chunk, not upon completion of each subchunk
	if len(conformer_list) > max_ligands_to_keep:
		conformer_heap = conformer_list

		print("Heapifying list after end of current file")
		heapq.heapify(conformer_heap)

		#use the nlargest command to get max number back as a list to use for the next file
		conformer_list = heapq.nlargest(max_ligands_to_keep,conformer_heap)

		#pop off until the length is at the max
		#while len(conformer_heap) > max_ligands_to_keep:
		#	heapq.heappop()

#write the best ligands to a file with accession information to a new file

write_file = open(str(file_prefix) + "_best_" + str(max_ligands_to_keep) + "_chunks_" + str(min_chunk_str) + "_" + str(max_chunk_str) + ".txt","w")

#print(list(conformer_list[i]))
#print("-------------")
#print()

for i in conformer_list:
	write_file.write(str(i[0]) + "," + str(i[1]) + "," + str(i[2]) + "," + str(i[3]) + "\n")

write_file.close()