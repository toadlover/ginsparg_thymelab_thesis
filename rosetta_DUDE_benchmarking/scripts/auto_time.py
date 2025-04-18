import os

os.mkdir("time_out")
os.system('find . -name "*.out" -exec cp "{}" time_out  \;')


os.chdir("time_out")
w = open("time_total.txt", "w")
for f1 in os.listdir("."):
    if (not f1.startswith("time")) & (not f1.startswith(".")):
        with open(f1, 'rb') as f:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
            last_line = f.readline().decode()
            w.write(f1 + "\n")
            w.write(last_line.split()[-1] + "\n\n\n")
            f.close()
w.close()


