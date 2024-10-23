import os,sys

#make dictionaries to hold the contents of the csv files
conf_counter_dict = {}

elements_counter_dict = {}

mw_no_h_dict = {}

#make dictionaries that get average values for MW and conf num for each chunk (help observe trends across the library)
average_conf_per_chunk_dict = {}

average_mw_no_h_per_chunk_dict = {}

#iterate over the entire conformer library
for i in range(0,53085):
#for i in range(0,10):
	#build the accession string
	chunk_str = str(i)
	while len(chunk_str) < 5:
		chunk_str = "0" + chunk_str

	#pull all csvs
	#os.system("s3cmd get s3://ariosg/ligand_library/" + chunk_str + "/for_s3/*.csv")
	os.system("s3cmd get s3://ariosg/ligand_library/" + chunk_str + "/for_s3/elements_encountered.csv s3://ariosg/ligand_library/" + chunk_str + "/for_s3/molecular_weight_no_h.csv")

	#initially seed average per chunk as 0, in case there in an issue with the chunk
	average_conf_per_chunk_dict[chunk_str] = 0
	average_mw_no_h_per_chunk_dict[chunk_str] = 0

	#read the csvs and add their contents to the counts, making new keys as necessary

	#conf count
	try:
		readfile = open("conformer_counter.csv", "r")

		#variables to hold the number of ligands observed and total number of conformers
		#average conformers in chunk derived as total conf num/number of ligands observed
		total_lig = 0
		total_conf = 0

		for line in readfile.readlines():
			conf_num = int(line.split(",")[0])
			count = int(line.split(",")[1].strip())

			#apply values to conf counter dictionary
			#add to count if key already exists, otherwise define the new key with value as current count if encountered for first time
			if conf_num in conf_counter_dict.keys():
				conf_counter_dict[conf_num] = conf_counter_dict[conf_num] + count
			else:
				conf_counter_dict[conf_num] = count

			#apply total lig and total conf to variables
			total_lig = total_lig + count
			total_conf = total_conf + (count * conf_num)

		#derive the average number of confs per ligand
		average_confs_per_ligand = total_conf/total_lig

		average_conf_per_chunk_dict[chunk_str] = average_confs_per_ligand
		readfile.close()
	except:
		print("conformer count for " + chunk_str + " is missing")

	#mw count
	try:
		readfile = open("molecular_weight_no_h.csv", "r")

		#variables to hold the number of mass bin and ligand count
		#average mass in chunk derived as mass * # ligands/total number of ligands observed
		total_lig = 0
		total_mass = 0

		for line in readfile.readlines():
			mass = int(line.split(",")[0])
			count = int(line.split(",")[1].strip())

			#apply values to conf counter dictionary
			#add to count if key already exists, otherwise define the new key with value as current count if encountered for first time
			if mass in mw_no_h_dict.keys():
				mw_no_h_dict[mass] = mw_no_h_dict[mass] + count
			else:
				mw_no_h_dict[mass] = count

			#apply total lig and total conf to variables
			total_lig = total_lig + count
			total_mass = total_mass + (count * mass)

		#derive the average number of confs per ligand
		average_mass_per_ligand = total_mass/total_lig

		average_mw_no_h_per_chunk_dict[chunk_str] = average_mass_per_ligand
		readfile.close()
	except:
		print("molecular weight for " + chunk_str + " is missing")

	#element count
	try:
		readfile = open("elements_encountered.csv", "r")

		for line in readfile.readlines():
			element = str(line.split(",")[0])
			count = int(line.split(",")[1].strip())

			#apply values to conf counter dictionary
			#add to count if key already exists, otherwise define the new key with value as current count if encountered for first time
			if element in elements_counter_dict.keys():
				elements_counter_dict[element] = elements_counter_dict[element] + count
			else:
				elements_counter_dict[element] = count
		readfile.close()
	except:
		print("elements count for " + chunk_str + " is missing")

	#remove csv files
	os.system("rm *csv")

#we now have all data from all csvs, print out to new csvs with data on whole library

#conformer count
writefile = open("all_conformer_count.csv", "w")

writefile.write("Conformers Generated, Ligands\n")

for confs in conf_counter_dict.keys():
	writefile.write(str(confs) + "," + str(conf_counter_dict[confs]) + "\n")

writefile.close()

#average conf per chunk
writefile = open("average_conformers_per_chunk.csv", "w")

writefile.write("Chunk, Average Conformers Generated\n")

for chunk in average_conf_per_chunk_dict.keys():
	writefile.write(str(chunk) + "," + str(average_conf_per_chunk_dict[chunk]) + "\n")

writefile.close()

#molecular weight
writefile = open("all_mw_count.csv", "w")

writefile.write("Molecular Weight (amu), Ligands\n")

for weight in mw_no_h_dict.keys():
	writefile.write(str(weight) + "," + str(mw_no_h_dict[weight]) + "\n")

writefile.close()

#average weight per chunk
writefile = open("average_mw_per_chunk.csv", "w")

writefile.write("Chunk, Average Molecular Weight (amu)\n")

for chunk in average_mw_no_h_per_chunk_dict.keys():
	writefile.write(str(chunk) + "," + str(average_mw_no_h_per_chunk_dict[chunk]) + "\n")

writefile.close()

#elements encountered
writefile = open("all_elements_count.csv", "w")

writefile.write("Element, Total Count\n")

for element in elements_counter_dict.keys():
	writefile.write(str(element) + "," + str(elements_counter_dict[element]) + "\n")

writefile.close()
