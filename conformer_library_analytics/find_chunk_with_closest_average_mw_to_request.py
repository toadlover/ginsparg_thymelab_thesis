#The purpose of this script is to input a target molecular weight, and then see what chunk in the conformer_library_analytics/average_mw_per_chunk.csv has the closest mw to it to be used in the benchmarking against other papers study

#the mw of interest is taken in as a command line argument, i.e.
#python find_chunk_with_closest_average_mw_to_request.py 297.113

#import
import os,sys

#take in the mw of interest
mw_of_interest = float(sys.argv[1])

#variables to hold the chunk with the closest average MW, the corresponding MW, and difference from input
closest_chunk_name = ""
closest_chunk_mw = ""
closest_chunk_diff = ""

#open and read through the average mw per chunk file
average_chunk_file = open(os.path.dirname(os.path.abspath(__file__)) + "/" + average_mw_per_chunk.csv, "w")

for line in average_chunk_file.readlines():
	#skip first line
	if "Average Molecular Weight" in line:
		continue

	#strip the newline off the line
	stripped_line = line.strip()

	#extract the chunk name and average mw
	cur_chunk_name = stripped_line.split(",")[0]
	cur_chunk_mw = float(stripped_line.split(",")[0])

	#get the absolute difference of the chunk mw and input mw
	cur_chunk_diff = abs(cur_chunk_mw - mw_of_interest)

	#assign if on first chunk
	if closest_chunk_diff == "":
		closest_chunk_diff = cur_chunk_diff
		closest_chunk_name = cur_chunk_name
		closest_chunk_mw = cur_chunk_mw
		continue

	#check if the difference is less, keep if less
	if cur_chunk_diff < closest_chunk_diff:
		closest_chunk_diff = cur_chunk_diff
		closest_chunk_name = cur_chunk_name
		closest_chunk_mw = cur_chunk_mw

#print the closest chunk
print("Closest chunk name: " + closest_chunk_name)
print("Closest chunk average MW: " + str(closest_chunk_mw))
print("Closest chunk diff from input: " + str(closest_chunk_diff))

