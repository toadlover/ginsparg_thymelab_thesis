import os,sys

for i in range(0,53084):

	#make a 5 digit string out of i
	i_str = str(i)

	while len(i_str) < 5:
		i_str = "0" + i_str

	print(i_str)

	os.system("mkdir output/" + i_str)
	os.system("mkdir error/" + i_str)
	os.system("mkdir log/" + i_str)
