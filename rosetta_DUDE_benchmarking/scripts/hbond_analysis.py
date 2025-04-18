import os
import statistics
#20 for attr, 15 for Hbond
os.chdir("rmsd_out")
w = open("hbond_analysis.txt", "w")
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
g = {"srbb": [], "lrbb": [], "bbsc": [], "sc": []}
l = {"srbb": [], "lrbb": [], "bbsc": [], "sc": []}
for k, v in m.items():
    if len(k.split("_")) < 15:
        continue
    v = v.replace("\n", "")
    k = k.replace(".pdb", "")
    if float(v) > 5:
        k = k.split("_")
        g.get("srbb").append(float(k[20]))
        g.get("lrbb").append(float(k[22]))
        g.get("bbsc").append(float(k[24]))
        g.get("sc").append(float(k[26]))
    else:
        k = k.split("_")
        l.get("srbb").append(float(k[20]))
        l.get("lrbb").append(float(k[22]))
        l.get("bbsc").append(float(k[24]))
        l.get("sc").append(float(k[26]))

w.write("len, mean, stddev, min, max\n")
w.write("> 5A\n")
w.write("\nsrbb\n")
w.write(str([len(g.get("srbb")), statistics.mean(g.get("srbb")), statistics.stdev(g.get("srbb")), min(g.get("srbb")), max(g.get("srbb"))]))
w.write("\nlrbb\n")
w.write(str([len(g.get("lrbb")), statistics.mean(g.get("lrbb")), statistics.stdev(g.get("lrbb")), min(g.get("lrbb")), max(g.get("lrbb"))]))
w.write("\nbbsc\n")
w.write(str([len(g.get("bbsc")), statistics.mean(g.get("bbsc")), statistics.stdev(g.get("bbsc")), min(g.get("bbsc")), max(g.get("bbsc"))]))
w.write("\nsc\n")
w.write(str([len(g.get("sc")), statistics.mean(g.get("sc")), statistics.stdev(g.get("sc")), min(g.get("sc")), max(g.get("sc"))]))
w.write("\n\n")
w.write("<= 5A\n")
w.write("\nsrbb\n")
w.write(str([len(l.get("srbb")), statistics.mean(l.get("srbb")), statistics.stdev(l.get("srbb")), min(l.get("srbb")), max(l.get("srbb"))]))
w.write("\nlrbb\n")
w.write(str([len(l.get("lrbb")), statistics.mean(l.get("lrbb")), statistics.stdev(l.get("lrbb")), min(l.get("lrbb")), max(l.get("lrbb"))]))
w.write("\nbbsc\n")
w.write(str([len(l.get("bbsc")), statistics.mean(l.get("bbsc")), statistics.stdev(l.get("bbsc")), min(l.get("bbsc")), max(l.get("bbsc"))]))
w.write("\nsc\n")
w.write(str([len(l.get("sc")), statistics.mean(l.get("sc")), statistics.stdev(l.get("sc")), min(l.get("sc")), max(l.get("sc"))]))
w.close()



