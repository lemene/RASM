
import os
import sys


quast_file = sys.argv[1]
bed = sys.argv[2]

with open(quast_file, "r") as f, open(bed, "w") as fo:
    for line in f:
        if line.startswith("chr") or line.startswith("ctg"):
            chr = line.strip()
            continue
        # 
        ls = line.strip().split(" ")
        # print(ls)
        # break
        for i in range(len(ls)):
            if ls[i] == "between":
                # print(ls[i+1], ls[i+2], ls[i+4], ls[i+5])
                reg1_ls = [int(ls[i+1]), int(ls[i+2])]
                reg2_ls = [int(ls[i+4]), int(ls[i+5])]
                fo.write("{}\t{}\t{}\n".format(chr, min(reg1_ls), max(reg1_ls)))
                print()
                break
            else:
                continue
