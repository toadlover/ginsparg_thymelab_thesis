#basic script to break up bigger files into more bite-sized chunks for parallel operations
#this will break the file up into 100 line segments, completely agnostic of the data
#each resulting file will be moved into a numbered directory, starting at zero
#the directory number will also be added to the front of the smaller file name

import os,sys

#read in the file of interest
read_file = open(sys.argv[1], "r")

#declare a line counter and file counter
line_counter = 0
file_counter = 0

#make a directory for file 0
os.system("mkdir " + str(file_counter))

#open write stream for first file
write_file = open(str(file_counter) + "/" + str(file_counter) + "_" + sys.argv[1], "w")

for line in read_file.readlines():
    write_file.write_file(line)

    #increment the line counter
    line_counter = line_counter + 1

    #if the line counter is divisible by 100, move to the next file
    if line_counter % 100 == 0:
        file_counter = file_counter + 1
        os.system("mkdir " + str(file_counter))

        #open write stream for first file
        write_file = open(str(file_counter) + "/" + str(file_counter) + "_" + sys.argv[1], "w")
                