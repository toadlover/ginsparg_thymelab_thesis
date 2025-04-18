import os
import statistics
w = open("stats.txt", "w")

root = "/data/user/abgvg9/benchmark/files/all"
initial_lig = 0
for item in os.listdir(root):
    if os.path.isdir(os.path.join(root, item)):
        initial_lig += 1

path, dirs, files = next(os.walk("/data/user/abgvg9/benchmark/motif_collection/files/neighbor_out"))
neighbor = len(files)

path, dirs, files = next(os.walk("disc_out"))
num_aa = len(dirs)

path, dirs, files = next(os.walk("rmsd_out"))
rmsd = len(files)

l = []
count = 0
total_c = 0
o = open("rmsd_out/rmsd_total.txt", "r")
for i in o.readlines():
    count += 1
    if count == 2:
        l.append(float(i))
    if count == 4:
        count = 0
        total_c += 1
t = []
time_out = 0
count = 0
o = open("time_out/time_total.txt", "r")
for i in o.readlines():
    count += 1
    if count == 2:
        if i.replace('.','').replace("\n", "").isnumeric():
            t.append(float(i))
        else:
            time_out += 1
    if count == 4:
        count = 0
t = [x / 60.0 for x in t]

w.write("Initial Receptors \n" + str(initial_lig))
w.write("\nNumber of Params Generated\n" + str(neighbor))
w.write("\nNumber of AA found\n" + str(num_aa))
w.write("\nNumber of Residues with placements found\n" + str(rmsd))
w.write("\nNumber of Receptors with placements found\n" + str(total_c))
w.write("\n\nPlacements within 1A\n" + str(sum(x <= 1 for x in l)))
w.write("\nPlacements between 1A and 2A\n" + str(sum(((x > 1) & (x <= 2)) for x in l)))
w.write("\nPlacements between 2A and 5A\n" + str(sum(((x > 2) & (x <= 5)) for x in l)))
w.write("\nPlacements greater than 5A\n" + str(sum(x >= 5 for x in l)))
w.write("\n\nAvg Time: " + str(statistics.mean(t)))
w.write("\nMax Time: " + str(max(t)))
w.write("\nTime within 5: " + str(sum(x <= 5 for x in t)))
w.write("\nTime between 5 and 10: " + str(sum(((x > 5) & (x <= 10)) for x in t)))
w.write("\nTime between 10 and 20: " + str(sum(((x > 10) & (x <= 20)) for x in t)))
w.write("\nTime greater than 20: " + str(sum(x > 20 for x in t)))
w.write("\nTimeouts: " + str(time_out))

w.close()



