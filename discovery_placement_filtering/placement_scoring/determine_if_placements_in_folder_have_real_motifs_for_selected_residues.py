#the purpose of this script is to look at a selected folder, read in a comma-separated list of residue indices (no spaces) that are wanted to check for if there are real motifs, and then print a list of the pdbs with the real motifs
#residue indices should be as they appear in the the Rosetta real motif code
#i.e. for the following line:
#Placement motifs: Real motif check - 7l1u_receptor_only_ResPos_86_ResID_GLN_Trio32_PV-005799603862_11_motif_7CTQ_FAD.pdb_GLN195_FAD_Packing_score:-1.381385_Hbond_score:0.000000_PRO131_PV-_Packing_score:-1.796332_Hbond_score:0.000000 remark: 6ZKB_PC1_Packing_score:-1.569930_Hbond_score:0.000000, distance: 0.709431, angle: 0.245363

#We are only concerned with the index value at PRO131, which corresponds to the motif found and taken off of the placement


#example call: python determine_if_placements_in_folder_have_real_motifs_for_selected_residues.py /data/project/thymelab/agonist_gln_results_pulldown/0/placements 134,118

#imports
import os,sys


#take in command line arguments
working_location = sys.argv[1]
residue_indices = sys.argv[2]

#interpret the residue list into indices to be iterated over
residue_list = residue_indices.split(",")



#create a checklist list that turns each residue index into a tuple of the index (as a string) and boolean to confirm if the residue has a real motif; defauly is false in case no interaction was even found, and will flip true if a real motif exists
residue_checklist = []
for item in residue_list:
	residue_checklist.append([item,False])

#make a copy of the list that can be used to quickly revert the working checklist
#residue_checklist_copy = residue_checklist

print("starting checklists")
print("original: ", residue_checklist)
#print("copy: ", residue_checklist_copy)

#build a string for the write file that has all requested indices listed in the name (used to differentiate between runs on the same folder if needed)
write_file_name = "real_motif_files"

for item in residue_list:
	write_file_name = write_file_name + "_" + item

write_file_name = write_file_name + ".txt"

#make a file to write the files with paths that contain real motifs for the requested residues
write_file = open(write_file_name, "w")

#move to the working location and then iterate over all files within it, only looking at pdb placement files that are at the working level
os.chdir(working_location)

for r,d,f in os.walk(working_location):
	for file in f:
		#confirm the file is what we want
		if file.endswith(".pdb") and r == working_location:

			print("on file: " + file)

			#wipe the checklist
			#residue_checklist = residue_checklist_copy

			for i in range(len(residue_checklist)):
				residue_checklist[i][1] = False

			print("Pre:")
			print("original: ", residue_checklist)
			#print("copy: ", residue_checklist_copy)

			#begin to read the file
			read_file = open(r + "/" + file,"r")

			for line in read_file.readlines():
				#we are only concerned with lines that contain "Real motif check"
				if "Real motif check" not in line:
					continue

				#if "remark:" isn't in the line, then the motif is fake too, so continue
				if "remark:" not in line:
					continue

				#if we have a line we want, split, and the motif info block is on the 7th "word"
				motif_str = line.split()[6]

				#extract the residue, which is 6th from the back when splitting by underscores
				residue = motif_str.split("_")[len(motif_str.split("_")) - 6]

				#extract the index number, cut off the first 3 characters from residue, which is the residue 3 letter code
				current_index = residue[3:]

				print(current_index)

				#iterate over the checklist and see if the index exists, and label true if so
				for i in range(len(residue_checklist)):
					if residue_checklist[i][0] == current_index:
						residue_checklist[i][1] = True

			print("Post:")
			print("original: ", residue_checklist)
			#print("copy: ", residue_checklist_copy)

			#bool to determine whether to keep the placement
			keep_placement = True

			#we have looked over the file, determine if if has real motifs for the desired indices
			for item in residue_checklist:

				if item[1] == False:
					keep_placement = False
					break

			#write the placement if we want to keep it
			if keep_placement:
				write_file.write(r + "/" + file + "\n")
				print("File has all desired real motifs!")

