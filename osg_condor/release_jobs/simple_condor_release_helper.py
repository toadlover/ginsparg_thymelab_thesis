#this script helps to release jobs on a regular basis as long as there are jobs in the condor queue (condor_q)
#the user inputs a time (in seconds) for how long the script checks condor_q to attempt to release

#run as "python3 simple_condor_release_helper.py 60 &"
#change the "60" to whatever refresh rate you want in seconds (ideally do not make it lower than like 10 seconds to minimize obtrusiveness on condor)

import os,sys

#rate to check condor_q (in seconds; ideally 5+ minutes, 300+ seconds)
refresh_rate = str(sys.argv[1])

#bool to indicate whether to keep checking
keep_checking = True

#get your user id for use in running condor_release and knowing when to end
os.system("whoami > whoami.txt")

#grab your id
read_file = open("whoami.txt", "r")

my_id = ""

for line in read_file.readlines():
	my_id = line.strip()

while keep_checking:
	#set keep checkign to False, we will set it to true if there are jobs that are still running
	keep_checking = False

	#pull current condor_q to read
	os.system("condor_q > condor_q_state.txt")

	#read condor_q_state
	read_file = open("condor_q_state.txt", "r")

	for line in read_file.readlines():
		#if "HOLD" is in any line (would only be in one of the condor_q header lines), run condor_release on my id (abgvg9)
		if "HOLD" in line:
			os.system("condor_release " + my_id)

		#if my id is in any lines (would only be in lines that correspond to a job), set keep_checking to True so that we can keep checking until all jobs are done
		if my_id in line:
			keep_checking = True

	#end by printing the date and condor_q state and then sleep for the refresh rate
	os.system("condor_q")
	os.system("date")
	os.system("echo SLEEPING FOR " + refresh_rate + " seconds")
	os.system("sleep " + refresh_rate)



#end
print("All jobs are done, bye!")

os.system("rm condor_q_state.txt whoami.txt")
