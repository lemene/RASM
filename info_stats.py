from multiprocessing import Pool
import os
import math
import pysam
from create_consensus_by_bed import fasta_parser
from collections import defaultdict


def cal_Ne(work_dir, ctg, ctg_len, min_dp, mm_rate):
    Ne, Na = 0, 0
    
    region = ctg + ":" + str(0) + "-" + str(ctg_len)
    pileup_file = os.path.join(work_dir, "pileup", region + ".pileup.txt")
    mis_bed = os.path.join(work_dir, "filtered2", region + ".bed")
    mis_ls = []
    with open(mis_bed, "r") as f:
        for line in f:
            fields = line.strip().split()
            mis_ls.append([fields[0], int(fields[1]), int(fields[2])])        
    if len(mis_ls) > 0:
        mis = mis_ls[0]
    else:            
        mis = [ctg, -1, -1]
    mis_length = 0
    mis_length += mis[2] - mis[1]
    mis_idx = 0
    
    with open(pileup_file, "r") as f:
        for line in f:
            items = line.strip().split("\t")
            pos = int(items[1])
            if pos >= mis[1] and pos <= mis[2]: 
                continue 

            if pos > mis[2]:  # over this mis
                if mis_idx + 1 < len(mis_ls):
                    mis_idx += 1
                    mis = mis_ls[mis_idx]
                    mis_length += mis[2] - mis[1]

            depth = int(items[3])
            match = items[4].count('.') + items[4].count(',')
            if depth >= min_dp:
                Na += 1
                if match / depth < mm_rate:
                    Ne += 1
    return Ne, Na

def get_qv(work_dir, ctgs, threads):

    simple_file = os.path.join(work_dir, "simple.ststs")
    f = open(simple_file, "a+")
    min_dp = 5
    mm_rate = 0.7   # 根据读数
    # 
    results = []
    pool = Pool(processes=threads)
    for ctg, ctglen in ctgs:
        results.append(pool.apply_async(cal_Ne, args=(work_dir, ctg, ctglen, min_dp, mm_rate)))

    pool.close() 
    pool.join()

    Ne = 0
    Na = 0
    for res in results:
        ctg_Ne, ctg_Na = res.get()
        Ne += ctg_Ne
        Na += ctg_Na
    qv = -10 * math.log10(Ne/Na)

    f.write("QV: {}\n".format(qv))
    f.close()