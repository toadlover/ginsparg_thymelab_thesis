#the purpose of this job is to write a queue file for condor usage that covers all subchunks, listing their chunk and superchunk
#this is for a broad discovery approach being developed in May 2025 where discovery is run on select ligands all in given subchunks, and so you can access the whole subchunk
#this writes a csv file that lists the superchunk,chunk,subchunk

out_file = open("all_subchunks_list.csv", "w")

#write the subchunk data based on loop iterations
for i in range(53085):
	#derive the superchunk as the floor of i/100
	superchunk = str(int(i/100))

	#derive the chunk as i
	chunk = str(i)

	#append leading zeroes for low number chunks
	while len(chunk) < 5:
		chunk = "0" + chunk

	for j in range(10):


		#derive the subchunk
		subchunk = str(j)

		#write the line to the csv for the chunk
		out_file.write(superchunk + "," + chunk + "," + subchunk)