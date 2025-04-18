import os
import statistics
#20 for attr, 15 for Hbond
os.chdir("rmsd_out")
w = open("faatr_analysis.txt", "w")
m = {}
for f1 in os.listdir("."):
    if (not f1.startswith("rmsd")) & (not f1.startswith(".")):
        o = open(f1, 'r')
        b = False
        l = []
        for line in o.readlines():

            if b:
                l.append(line.replace("\n", ""))
                m[l[0]] = l[1]
                l = []

                b = False
            elif line.startswith("no"):
                b = True
                l.append(line)
        o.close()
g = []
l = []
for k, v in m.items():
    if len(k.split("_")) < 15:
        continue
    v = v.replace("\n", "")
    if float(v) > 5:
        temp = k.split("_")[15]
        g.append(float(temp))
    else:
        temp = k.split("_")[15]
        l.append(float(temp))

w.write("len, mean, stddev, min, max for > 5A\n")
w.write(str([len(g), statistics.mean(g), statistics.stdev(g), min(g), max(g)]))
w.write("\n\n")
w.write("len, mean, stddev, min, max for <= 5A\n")
w.write(str([len(l), statistics.mean(l), statistics.stdev(l), min(l), max(l)]))

w.close()



