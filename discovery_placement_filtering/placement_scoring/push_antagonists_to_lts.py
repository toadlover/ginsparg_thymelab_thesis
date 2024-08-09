import os,sys

for r,d,f in os.walk(os.getcwd()):
	for file in f:
		if file == "test_params.tar.gz":

			#get the directory
			dir_num = str(r.split("/")[len(r.split("/")) - 1])

			print("s3cmd put " + r + "/test_params.tar.gz s3://ariosg/antagonists_12M_chunk_sorted/" + dir_num + "/ --no-progress")