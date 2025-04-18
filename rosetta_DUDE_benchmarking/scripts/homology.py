import os
import re

o = open("/data/user/abgvg9/benchmark/files/fasta/fasta_out.log", "r")
m_f = {}
count = 0
cur_p = ""
for i in o.readlines():
    i = i.strip()
    if (not i) & (count == 8):
        count = 0
    elif count == 8:
        prot = re.search('(.*)_', i).group(1)
        m_f.get(cur_p).append(prot)
    elif count != 0:
        count += 1
    elif i.startswith("Query="):
        prot = re.search('Query= (.*)_', i).group(1)
        m_f[prot] = []
        cur_p = prot
        count += 1
o.close()
m_d = {}
for root, dirs, files in os.walk("disc_out"):
    for file in files:
        if (not file.endswith("lig.pdb")) & (file.startswith("noclash")):
            p = os.path.join(root, file)
            s = p.split("_")
            pdb = s[11]
            prot = re.search('/(.*)', s[1]).group(1)
            if not prot in m_d:
                m_d[prot] = []
            else:
                m_d.get(prot).append(pdb)
m_r = {}
r = open("rmsd_out/rmsd_total.txt", "r")
count = 0
for i in r.readlines():
    if i.startswith("noclash"):
        s = i.split("_")
        pdb = s[7]
        prot = s[1]
        m_r[prot] = pdb

o = open("homology_out.txt", "w")
for key1, value1 in m_d.items():
    l_d = value1
    l_f = m_f.get(key1)
    sim = []
    for i1 in l_d:
        for i2 in l_f:
            if (i1.lower() == i2.lower()) & (i1 not in sim):
                sim.append(i1)
    o.write(key1 + ": \n")
    for temp in sim:
        o.write("\t" + temp)
        print(m_r.get(key1), key1)
        if m_r.get(key1).lower() == temp.lower():
            o.write(" <- min RMSD ")
    o.write("\n")

o.close()
