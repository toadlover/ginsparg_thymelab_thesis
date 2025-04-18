import os
import time
import re
os.system("source activate rosetta_tools_env")
from datetime import date, datetime

argm = {}
o = open("/data/user/abgvg9/benchmark/bench_params", "r")
for line in o.readlines():
        line = line.split()
        argm[line[0]] = line[1]
o.close()

def checkJobs():
    os.system("squeue -u " + user + " > squeue_length.txt")
    s = open("squeue_length.txt")
    c = 0
    for line in s.readlines():
        if not line.split()[2].startswith("Inter"):
            c += 1
    os.system("rm squeue_length.txt")

    return c > 2



os.chdir("/data/user/abgvg9/benchmark")

count = 0
user = ""
td = str(datetime.today().strftime('%m-%d-%Y'))

os.system("whoami > user.txt")
o = open("user.txt")
user = o.readlines()[0].strip()
o.close()
os.system("rm user.txt")
for f in os.listdir("runs"):
    if f.startswith("run"):
        ctemp = re.search("run(.*)_", f).group(1)
        print(ctemp)
        if int(ctemp) > int(count):
            count = int(ctemp)
count += 1

os.mkdir('runs/run' + str(count) + "_" + td)
dir_loc = 'runs/run' + str(count) + "_" + td + "/"

o = open(dir_loc + "run" +str(count) + "_info.txt", "w")


o.write("User: " + user + "\n")
o.write("Arguments: " + str(argm) + "\n")

os.chdir(dir_loc)

print("run" + str(count) + "\n")
print("Info Created\n")
o.write("Info Created\n\n")

print("Placing Ligands\n")
o.write("Placing Ligands\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")
st = str(datetime.now())
os.system("python3 /data/user/abgvg9/benchmark/scripts/auto_place.py " +
          argm["-l"] + " " + argm["-b"] + " " + argm["-r"] + " " + argm["-a"] + " " + argm["-s"] + " " + argm["-m"]  + " " + argm["-d"])

time.sleep(15)

while checkJobs():
    time.sleep(15)

et = str(datetime.now())
print("Ligands Placed\n\n")
o.write("Ligands Placed\n\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")
print("Calculating RMSD\n")
o.write("Calculating RMSD\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")
os.system("python3 /data/user/abgvg9/benchmark/scripts/auto_rmsd.py")
print("RMSD Done\n\n")
o.write("RMSD Done\n\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")

print("Performing RMSD Analysis\n")
o.write("Performing RMSD Analysis\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")
os.system("python3 /data/user/abgvg9/benchmark/scripts/auto_rmsd_total.py")
print("RMSD Calculations Done\n\n")
o.write("RMSD Calculations Done\n\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")

print("Performing Time Analysis\n")
o.write("Performing Time Analysis\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")
os.system("python3 /data/user/abgvg9/benchmark/scripts/auto_time.py")
print("Time Calculations Done\n\n")
o.write("Time Calculations Done\n\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")

print("Performing Statistical Analysis\n")
o.write("Performing Statistical Analysis\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")
os.system("python3 /data/user/abgvg9/benchmark/scripts/stats.py")
print("Statistical Calculations Done\n\n")
o.write("Statistical Calculations Done\n\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")

print("Performing Homology Analysi\n")
o.write("Performing Homology Analysis\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")
os.system("python3 /data/user/abgvg9/benchmark/scripts/homology.py")
print("Homology Calculations Done\n\n")
o.write("Homology Calculations Done\n\n")
o.close()
o = open("run" +str(count) + "_info.txt", "a")

o.write("Start Time: " + st + "\n")
print("Start Time: " + st + "\n")
print("End Time: " + et)
o.write("End Time: " + et)
o.close()




