import sys

fout = sys.argv[1]
i = 0
with open(fout, "w") as f:
    f.write("{}\t{}\n".format(i, ))