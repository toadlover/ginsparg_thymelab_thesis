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
						if control_group in dire and experimental_group in dire and r == expt[0]:
							os.system("cp " + r + "/" + dire + "/ribgraph_mean_" + section + "_" + metric + "_" + my_bin + ".png .")
							os.system("cp " + r + "/" + dire + "/boxgraph_ribgraph_mean_" + section + "_" + metric + "_" + my_bin + ".png_data.csv .")

							#modify the csv file so that the control group column is always listed first (since it is possible that it could be swapped)
							#also change column labels so that they are exact matches of the experimental and control groups

							#first, determine if the control group is listed first or not
							#determine via seeing if the control group is before or after "_vs_" in the directory
							control_is_first = True
							if control_group in dire.split("_vs_")[1]:
								control_is_first = False

							#check if the file was suffessfully, and if not, continue so we do not break the script
							if os.path.exists("boxgraph_ribgraph_mean_" + section + "_" + metric + "_" + my_bin + ".png_data.csv") == False:
								continue

							#now, read through the file and modify it appropriately
							#write the modified file, and then copy it over the original
							orig_file = open("boxgraph_ribgraph_mean_" + section + "_" + metric + "_" + my_bin + ".png_data.csv", "r")
							write_file = open("temp.csv","w")

							#line counter
							line_counter = 0

							#read over the original and adjust the file
							for line in orig_file.readlines():
								#handle the header line
								if line_counter == 0:
									#write new header line
									write_file.write("\t" + control_group + "\t" + experimental_group + "\n")
								else:
									#other lines

									#break up the line by tabs
									split_line = line.split("\t")

									#reverse the data order if the control is listed second
									if control_is_first == False:
										flipped_line = [split_line[0], split_line[2], split_line[1]]
										split_line = flipped_line

									#print the potentially flipped line to the file
									write_file.write(split_line[0].strip() + "\t" + split_line[1].strip() + "\t" + split_line[2].strip() + "\n")

								#increment the line counter
								line_counter = line_counter + 1



							#close both files and then replace the original with the modded
							orig_file.close()
							write_file.close()
							os.system("mv temp.csv boxgraph_ribgraph_mean_" + section + "_" + metric + "_" + my_bin + ".png_data.csv")

	#at end, go up so we can do another directory or move onto the analysis
	os.chdir("..")

#Run analyses for each permutation of metric, section, and bins to make figures
#for each experiment metric
for metric in experiment_metrics:
	#for each time section:
	for section in time_sections:
		#for each bin
		for my_bin in bins:

			#this section was made with the help of chatgpt for figure wrangling
			# -----------------------------
			# Load and normalize data
			# -----------------------------
			base_dir = os.getcwd()
			pattern = os.path.join(base_dir, "*",  
			                       "boxgraph_ribgraph_mean_" + section + "_" + metric + "_" + my_bin + ".png_data.csv")
			files = glob.glob(pattern)
			print(f"Found {len(files)} data files.")

			all_data = []
			for file in files:
			    df = pd.read_csv(file, sep="\t")
			    
			    # Normalize relative to DMSO mean per experiment
			    dmso_mean = df[control_group].mean()
			    df[control_group] = df[control_group] / dmso_mean
			    df[experimental_group] = df[experimental_group] / dmso_mean
			    
			    # Add experiment metadata
#			    experiment_id = os.path.basename(os.path.dirname(os.path.dirname(file)))
			    experiment_id = os.path.basename(os.path.dirname(file))
			    
			    #if len(experiment_id.split("_")) > 1:
			    #    experiment_id = experiment_id.split("_")[1]
			    
			    df["experiment"] = experiment_id
			    
			    all_data.append(df[["experiment", control_group, experimental_group]])

			merged_df = pd.concat(all_data, ignore_index=True)

			# Sort merged_df alphabetically by experiment name
			merged_df = merged_df.sort_values(by="experiment").reset_index(drop=True)

			# -----------------------------
			# Convert to long format
			# -----------------------------
			plot_df = merged_df.melt(
			    id_vars="experiment", 
			    value_vars=[control_group, experimental_group],
			    var_name="treatment",
			    value_name="normalized_response"
			)
			plot_df["treatment"] = plot_df["treatment"].str.replace("_norm", "")

			# -----------------------------
			# Mann–Whitney U test per experiment
			# -----------------------------
			pvals = {}
			for exp in plot_df["experiment"].unique():
			    data_exp = plot_df[plot_df["experiment"] == exp]
			    dmso_vals = data_exp[data_exp["treatment"] == control_group]["normalized_response"].dropna()
			    drug_vals = data_exp[data_exp["treatment"] == experimental_group]["normalized_response"].dropna()
			    
			    if len(dmso_vals) >= 2 and len(drug_vals) >= 2:
			        _, p_val = mannwhitneyu(dmso_vals, drug_vals, alternative='two-sided')
			    else:
			        p_val = np.nan
			    pvals[exp] = p_val
			    
			# Print p-values
			for exp, p in pvals.items():
			    print(f"{exp}: p = {p}")

			# -----------------------------
			# Helper: p-value to stars
			# -----------------------------
			def p_to_stars(p):
			    if p < 0.001:
			        return "***"
			    elif p < 0.01:
			        return "**"
			    elif p < 0.05:
			        return "*"
			    else:
			        return "ns"

			# -----------------------------
			# Plot
			# -----------------------------
			plt.figure(figsize=(12,6))

			# Stripplot (points) with legend
			sns.stripplot(
			    data=plot_df,
			    x="experiment",
			    y="normalized_response",
			    hue="treatment",
			    dodge=True,
			    alpha=0.3,
			    palette={control_group: "black", experimental_group: "red"},
			    legend=False
			)

			# Barplot (mean ± SEM, outlined)
			sns.barplot(
			    data=plot_df,
			    x="experiment",
			    y="normalized_response",
			    hue="treatment",
			    estimator=np.mean,
			    errorbar=("se", 1),
			    dodge=True,
			    alpha=0,
			    fill=False,
			    palette={control_group: "black", experimental_group: "red"},
			    linewidth=2,
			    capsize=0.2,
			    errwidth=2,
			    legend=False
			)

			# Remove top/right spines
			sns.despine(trim=True)

			# Legend in top right
			#plt.legend(title="Treatment", loc="upper left")

			plt.xticks(rotation=45)
			plt.ylabel("Normalized " + metric + ": Binned - " + my_bin)
			plt.xlabel("")

			# Determine y positions for significance stars/brackets
			y_max_global = plot_df["normalized_response"].max()
			offset = 0.05  # vertical space above bars for bracket start

			for i, exp in enumerate(plot_df["experiment"].unique()):
			    p_val = pvals[exp]
			    stars = p_to_stars(p_val)
			    
			    if stars != "ns":
			        # Draw bracket line
			        bar_width = 0.2  # horizontal extent of bracket
			        x1 = i - bar_width
			        x2 = i + bar_width
			        y = y_max_global + 0.05
			        plt.plot([x1, x1, x2, x2], [y, y+0.02, y+0.02, y], color="black", lw=1.5)
			        
			        # Add stars above bracket
			        plt.text(i, y + 0.025, stars, ha='center', va='bottom', fontsize=14)

			# Move title below the plot
			plt.figtext(0.5, 0.01, section,
			            ha="center", fontsize=14)

			plt.tight_layout()
			#plt.show()

			#write the image
			plt.savefig(section + "_" + metric + "_" + my_bin + ".png", dpi=300, bbox_inches='tight', transparent=True)
			plt.close()