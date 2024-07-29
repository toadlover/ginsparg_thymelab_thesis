import os,sys

#call this script in a location where score csvs are
#looking for raw_scores.csv and weighted_scores.csv
#will append the working path to the front of the file names for easier file tracking (as long as the files aren't later moved from their locations)

#get the working location
location = os.getcwd()

#open the raw scores file and a temp file to write to
raw_file = open("raw_scores.csv", "r")
temp_file = open("temp.csv", "w")

for line in raw_file.readlines():
	#write header line as is
	if line.startswith("file,") and ",total" in line:
		temp_file.write(line)
	#write location to start of line and write to temp file
	else:
		temp_file.write(location + "/" + line)

#overwrite original raw file with filled temp
os.system("mv temp.csv raw_scores.csv")

#close ifle streams
raw_file.close()
temp_file.close()

#repeat process for weighted
#open the weighted scores file and a temp file to write to
weighted_file = open("weighted_scores.csv", "r")
temp_file = open("temp.csv", "w")

for line in weighted_file.readlines():
	#write header line as is
	if line.startswith("file,") and ",total" in line:
		temp_file.write(line)
	#write location to start of line and write to temp file
	else:
		temp_file.write(location + "/" + line)

#overwrite original raw file with filled temp
os.system("mv temp.csv weighted_scores.csv")

#close ifle streams
weighted_file.close()
temp_file.close()