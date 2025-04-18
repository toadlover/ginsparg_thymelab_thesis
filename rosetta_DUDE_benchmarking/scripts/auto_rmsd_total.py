import os

os.mkdir("rmsd_out")
os.system('find . -name "*rmsd_out.txt" -exec cp "{}" rmsd_out  \;')

v = []
os.chdir("rmsd_out")
w = open("rmsd_total.txt", "w")
for f1 in os.listdir("."):
    prot = f1.split("_")[0]
    if (not prot in v) & (not prot == "rmsd"):
        v.append(prot)
        l = []
        minR = 1000
        minL = ""
        for f2 in os.listdir("."):
            if f2.startswith(prot + "_"):
                o = open(f2, 'r')
                tminR = 0.0
                for line in o.readlines():
                    if line.startswith("BEST R"):
                        if (float(line.split(" ")[2]) < minR):
                            minR = float(line.split(" ")[2])
                            tminR = float(line.split(" ")[2])
                    elif line.startswith("BEST L") & (minR == tminR):
                        minL = line.split(" ")[2]
                o.close()
        w.write(prot + "\n" + str(minR) + '\n' + minL + "\n\n")

w.close()



