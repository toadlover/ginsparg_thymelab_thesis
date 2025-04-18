import os
import math
import re

def getRMSD(p1, p2):
    x = open(p1, "r")
    count = 0
    lx = []
    for line in x.readlines():
        if count > 6:
            line = line.strip()
            line = re.sub(' +',' ',line)
            line = line.split(" ")
            if len(line) < 9:
                break
            elif line[1].startswith("H"):
                continue
            lx.append((float(line[2]), float(line[3]), float(line[4])))
        count += 1

    ll = []
    count = 0
    x = open(p2, "r")
    for line in x.readlines():
        if line.startswith("HETATM"):
            line = line = re.sub(' +',' ',line)
            line = line.split(" ")
            if line[2].startswith("H"):
                continue
            ll.append((float(line[6]), float(line[7]), float(line[8])))


    sum = 0.0
    for i in range(len(ll)):
        sum += (lx[i][0] - ll[i][0]) ** 2
        sum += (lx[i][1] - ll[i][1]) ** 2
        sum += (lx[i][2] - ll[i][2]) ** 2

    sum /= len(lx)
    rmsd = math.sqrt(sum)
    return rmsd

os.chdir("disc_out")

counter = 0

dirlist = [ item for item in os.listdir(".") if os.path.isdir(os.path.join(".", item)) ]
for d in dirlist:
    
    counter = counter + 1
    
    print("RMSD",d,counter,len(dirlist))
    os.chdir(d)
    #moving these to be done in the place jobs to help speed things up (since this will be paralelized)
    os.system('''for file in *; do mv "$file" `echo $file | tr ' ' '_'` ; done''')
    os.system('''for file in *; do mv "$file" `echo $file | tr ':' '_'` ; done''')
    prot = d.split("_")[0]
    native = "/data/user/abgvg9/benchmark/files/all/" + prot + "/crystal_ligand.mol2"
    b = False
    m = {}
    for files in os.listdir("."):
        if "minipose" in files:
                continue
 #       if files.startswith("noclash") | files.startswith("debug") | files.endswith(".pdb"):
 #       if "_delta_" in files and files.endswith(".pdb"):
        #account for new and only naming conventions
        if (files.startswith(prot) and files.endswith(".pdb")) or ("_delta_" in files and files.endswith(".pdb")):
            b = True
            o = open(files ,'r')
            w = open(files.replace(".pdb", "") + '_lig.pdb','w')
            for line in o.readlines():
                if not (line.startswith('ATOM') | line.startswith("SSBOND") | line.startswith("HETNAM") ):
                    w.write(line)
            w.close()
            o.close()
            temp = files.replace(".pdb", "") + '_lig.pdb'
            m[files + "_lig"] = getRMSD(native, files.replace(".pdb", "") + '_lig.pdb')

    if b:
        w = open(d + "_rmsd_out.txt", "w")
        k = min(m, key = m.get)
        v = m[k]
        for x in m:
            w.write("\n" + str(x) + '\n' + str(m[x]) + '\n')
        w.write("\nBEST RMSD: " + str(v) + "\n" + "BEST LIG: " + str(k))

    os.chdir("..")


