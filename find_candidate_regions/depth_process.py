import os
import sys
import gzip
f_in = sys.argv[1]
f_out = sys.argv[2]

if f_in.endswith("gz"):
    f_reader = gzip.open(f_in, "rt")
else:
    f_reader = open(f_in, "r")
with open(f_out, "w") as f_writer:
    for line in f_reader:
        fields = line.strip().split()
        f_writer.write("{}\t{}\n".format(int(fields[1])//100, fields[3]))

f_reader.close()