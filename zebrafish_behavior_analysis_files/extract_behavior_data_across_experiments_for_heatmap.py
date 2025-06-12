#the purpose of this script is to look down all directories from a given location with a given substring for zebrafish behavior statistical analyses on given metrics 

#this script needs 3 inputs:

#1 - root location where all directories of interest lay

#2 - substring that all folders contain (used for selecting the right folders to look for data)
#i.e. "ari_data"

#3 - a PlotParameters file which contains all metrics that are desired to be looked at

#imports

import os,sys

#read in the inputs

root_location = sys.argv[1]

#have the root location end with a / if it doesn't already
if root_location.endswith("/") == False:
	root_location = root_location + "/"

shared_substring = sys.argv[2]
pp = sys.argv[3]

#first, collect all plot parameters of interest
pp_file = open(pp, "r")

#store parameters as a list
parameters = []
for line in pp_file.readlines():
	#sanity check, each parameter should be everything before the first :. No :, no parameter
	if ":" not in line:
		continue
	parameter = line.split(":")[0]
	parameters.append(parameter)

#now, run through each parameter and make a heatmap of data from all output data files (that have what we want)
for parameter in parameters:
	print(parameter)

	#list to hold drug data
	#each entry in the list will be a sublist that holds the full file name with path, day0night, day1night, day2night, day1day, day2day, night average, day average
	parameter_data = []

	#iterate over each out file
	for r,d,f in os.walk(root_location):
		for file in f:
			#make sure the shared substring is in the file root, and that the file starts with "linearmodel" and ends with ".out"
			if (shared_substring in r) and (file.startswith("linearmodel")) and (file.endswith(".out")):
				#we have a file to look at, concat a string of the root + file
				full_file = r + "/" + file

				#read the file and extract the terms for the parameter for the 3600 dpix measure
				#set variables with default values so we can make sure that we actually grab values from the file
				#if any value remains as X, then we failed to get values and should not post the data
				day0night = "X"
				day1night = "X"
				day2night = "X"
				day1day = "X"
				day2day = "X"
				night_avg = 0
				day_avg = 0

				#open the data file and try to collect the data
				data_file = open(full_file,"r")

				for line in data_file.readlines():
					#to consider a line, see if it starts with "anova:"
					if line.startswith("anova:") == False:
						continue

					#we are on an anova line, see if it has data we want
					#metric is line.split()[1]
					curr_metric = line.split()[1]

					#make sure it is dpix and 3600 (ends with either _3600 or _1over3600)
					if (curr_metric.endswith("_3600") or curr_metric.endswith("_1over3600")) and parameter in curr_metric:
						#extract the SSMD and assign it to the corresponding day, SSMD is the 3rd from last item in the line
						ssmd = line.split()[len(line.split())-3]

						#assign to appropriate value
						if "day0night" in curr_metric:
							day0night = ssmd
						if "day1night" in curr_metric:
							day1night = ssmd
						if "day2night" in curr_metric:
							day2night = ssmd
						if "day1day" in curr_metric:
							day1day = ssmd
						if "day2day" in curr_metric:
							day2day = ssmd

				#we now (hopefully) have the metric data. confirm we do, and if so, prepare to store it for writing to output csv
				#if any of the day/night values are still "X", continue and this was a bad file
				if day0night == "X" or day1night == "X" or day2night == "X" or day1day == "X" or day2day == "X":
					continue

				#keep the data

				#derive the night and day averages
				night_avg = (float(day0night) + float(day1night) + float(day2night)) / 3
				day_avg = (float(day1day) + float(day2day)) / 2

				drug_data = [full_file, day0night, day1night, day2night, day1day, day2day, night_avg, day_avg]
				parameter_data.append(drug_data)

	#we now have all the drug data for the parameter
	#write it to a csv file
	out_file = open(root_location + shared_substring + parameter + ".csv", "w")

	#write header line
	out_file.write("drug,day0night,day1night,day2night,day1day,day2day,night_avg,day_avg\n")

	#write drug lines
	for drug in parameter_data:
		out_file.write(str(drug[0]) + "," + str(drug[1]) + "," + str(drug[2]) + "," + str(drug[3]) + "," + str(drug[4]) + "," + str(drug[5]) + "," + str(drug[6]) + "," + str(drug[7]) + "\n")
