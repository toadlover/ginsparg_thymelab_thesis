#the purpose of this script it to take an input file that lists all locations of interest (as paths) to aggregate into figure bar images

#all included experiments must be the same type (i.e. knockout, wt, or heatshock) so that comparisons between the groups of interest can be made

#imports 
import sys
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

#make variables to hold the arguments to run this

#list of paths of all experiments for this study
#each entry will be a tuple of a path and a name to use to represent the experiment in the plot
experiment_paths = []

#bins and time_sections are used to select the csv files within the respective directories
#list of binning values to use (choice of 60,600,3600), must use 1, can use up to all 3
bins = []

#experiment time sections (i.e time_combo, day0night, etc.)
time_sections = []

#experiment metrics
experiment_metrics = []

#the 2 groups to compare (i.e. dmso-het vs drug-het for heat shock drug effect, dmso-wt vs dmso-het for heat shock confirmation, dmso-mut vs drug-mut to ko)
#this will help with searching for the right outputfulldata_... folder from each experiment
control_group = ""
experimental_group = ""

#read in the arguments file
#each line in the arguments file must start with the corresponding variable name followed by a colon for proper read-in
#for list variables, list an individual value on its own line, the experiment_paths tuple must be listed as comma-separated values for path,name
#argument file name can include the path, as is passed as the only command line argument for this
#order of arguments in the file does not matter
arg_file_name = sys.argv[1]
arg_file = open(arg_file_name,"r")

for line in arg_file.readlines():
	#experiment paths
	if line.startswith("experiment_paths:"):
		my_item = line.split("experiment_paths:")[1].strip()
		data_path = my_item.split(",")[0]
		data_name = my_item.split(",")[1]
		experiment_paths.append([data_path,data_name])
	#bins
	if line.startswith("bins:"):
		my_item = line.split("bins:")[1].strip()
		bins.append(my_item)
	#time sections
	if line.startswith("time_sections:"):
		my_item = line.split("time_sections:")[1].strip()
		time_sections.append(my_item)
	#experiment metrics
	if line.startswith("experiment_metrics:"):
		my_item = line.split("experiment_metrics:")[1].strip()
		experiment_metrics.append(my_item)
	#control group
	if line.startswith("control_group:"):
		my_item = line.split("control_group:")[1].strip()
		control_group = my_item
	#experiment group
	if line.startswith("experimental_group:"):
		my_item = line.split("experimental_group:")[1].strip()
		experimental_group = my_item

#make a working directory in the working location named after the arguments file that will copy all of the needed data into the folder so the figure can be made and an organized record kept
#name a directory based on what is before the first period in the argument file name, and add "_graphs" to the end
working_dir_name = arg_file_name.split(".")[0] + "_graphs"

#delete a previous directory and then make a new one, then move into it
os.system("rm -drf " + working_dir_name)
os.system("mkdir " + working_dir_name)
os.chdir(working_dir_name)

#debug print of the argument values
print("experiment_paths",experiment_paths)
print("bins",bins)
print("time_sections",time_sections)
print("experiment_metrics",experiment_metrics)
print("control_group",control_group)
print("experimental_group",experimental_group)

#now, for each experiment, make a sub-folder in the newly made folder and then work on copying all relevant data files (and images)
for expt in experiment_paths:
	#make a folder based on the name
	os.system("mkdir " + expt[1])
	#enter the folder
	os.chdir(expt[1])

	#for each experiment metric
	for metric in experiment_metrics:
		#for each time section:
		for section in time_sections:
			#for each bin
			for my_bin in bins:
				#copy the csv and line plot files that match the metric, section, and bin
				for r,d,f in os.walk(expt[0]):
					for dire in d:
						if control_group in dire and experimental_group in dire:
							os.system("cp " + r + "/" + dire + "/ribgraph_mean_" + section + "_" + metric + "_" + my_bin + ".png .")
							os.system("cp " + r + "/" + dire + "/boxgraph_ribgraph_mean_" + section + "_" + metric + "_" + my_bin + ".png_data.csv .")

	#at end, go up so we can do another directory or move onto the analysis