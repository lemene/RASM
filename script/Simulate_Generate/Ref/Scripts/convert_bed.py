# convert SV/mis bed from: ref
# To: ref+SV/mis

from collections import defaultdict
import os
import sys

bed_in = sys.argv[1]
bed_out = sys.argv[2]

class Bed:
    def __init__(self) -> None:
        pass
    @staticmethod
    def read_bed(bed):
        with open(bed, "r") as fin:
            dic = defaultdict(list)
            for line in fin:
                fileds = line.strip().split("\t")
                chr, start, end = fileds[:3]
                start, end = int(start), int(end)
                dic[chr].append((start, end))
            # sorted
            # chr_ls = list(dic.keys())
            # print(chr_ls)
            sorted_dic = defaultdict(list)
            for chr, reg_ls in dic.items():
                reg_ls.sort()
                sorted_dic[chr] = reg_ls
        return sorted_dic
        # pass
    @staticmethod
    def write_bed():
        pass


def read_bed(bed):
    with open(bed, "r") as fin:
        for line in fin:
            pass
    pass





