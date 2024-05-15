import os,sys

os.system(" condor_q -hold  | grep \"input files failure\" | awk \'{print $1}\' | grep 30675740 > missing_input_sdf.txt")
read_file = open("missing_input_sdf.txt")

for line in read_file.readlines():
	jobid = line.strip("\n")
	print(jobid)
	print("condor_rm " + str(jobid))
	os.system("condor_rm " + str(jobid))
