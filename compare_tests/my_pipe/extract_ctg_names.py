import sys
fai = sys.argv[1]
if fai.endswith("fai"):
    with open(fai, "r") as f:
        ctg_ls = []
        for line in f:
            ctg = line.split("\t")[0]
            ctg_ls.append(ctg)
    print(','.join(ctg_ls)) # contig_1,contig_11,contig_12,contig_14,contig_15,contig_16,contig_17,contig_19,contig_2,contig_20,contig_21,contig_23,contig_27,contig_28,contig_29,contig_3,contig_32,contig_35,contig_41,contig_43,contig_45,contig_52,contig_53,contig_54,contig_55,contig_58,contig_60,contig_63,contig_64,contig_66,contig_7,contig_8,contig_9
else:
    print("input fai file!!!")

