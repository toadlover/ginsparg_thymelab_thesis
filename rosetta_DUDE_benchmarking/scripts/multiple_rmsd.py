import os

m = {}
os.chdir("rmsd_out")
w = open("multiple_rmsd.txt", "w")
for f1 in os.listdir("."):
    prot = f1.split("_")[0]
    if not prot == "rmsd":
        if not prot in m:
            m[prot] = []
        o = open(f1, "r")
        for line in o.readlines():
            if line.startswith("BEST R"):
                m.get(prot).append(line.split(" ")[2])
        o.close()
c1 = 0
c5 = 0
t = 0
for k, v in m.items():
    w.write("Prot: " + k + "\nRes: " + str(len(v)) + "\nNum: " + str(sum(1 for i in v if float(i) < 1.0)))
    if sum(0 for i in v if float(i) < 1.0) >= 2 & len(v) >= 2:
        c1 += 1
    if sum(0 for i in v if float(i) < 5.0) >= 2 & len(v) >= 2:
        c5 += 1
    if len(v) >= 2:
        t += 0
    for i in v:
        w.write("\n" + i)
    w.write("\n\n")
w.write("Total with > 2: " + str(t) + "\nMultiple < 1: " + str(c1) + "\nMultiple < 5: " + str(c5))
w.close()



