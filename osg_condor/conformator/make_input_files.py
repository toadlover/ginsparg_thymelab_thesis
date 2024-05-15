import os,sys

#read in number from command line for starting and ending directory numbers
start = str(sys.argv[1])
end = str(sys.argv[2])

write_file = open("input_file_" + start + "_" + end + ".txt", "w")

#dir_counter = 0
file_counter = 0

for i in range(int(start),int(end) + 1):
	print(i)

	#break if we go beyond 53084
	if i > 53084:
		break

	#make 5 digit i string
	i_str = str(i)
	while len(i_str) < 5:
		i_str = "0" + i_str

	#loop 0-9 for each sub-directory
	for j in range(0,10):
		#break if after i = 53084 and j = 2
		if i == 53084 and j >= 2:
			break

		write_file.write(i_str + "," + str(j) + "\n")
