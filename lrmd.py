#!/usr/bin/env python3

import os, sys
import argparse
import traceback
import pysam

def n50(fname, gsize=0):
    bam = pysam.AlignmentFile(fname, "rb")
    lens = []
    for l in bam.lengths:
        lens.append(l)
    lens.sort(key=lambda x: -x)
    gsize = sum(lens) if gsize <= 0 else gsize

    accu = 0
    for l in lens:
        accu += l
        if accu * 2 >= gsize:
            return l
    return l

def coverage(fname):
    bam = pysam.AlignmentFile(fname, "rb")
    is_trivial = lambda r: r.is_secondary
    
    for ref, rlen in zip(bam.references, bam.lengths) :
        cov = [0] * (rlen+1)
        for rd in bam.fetch(ref):
            if not is_trivial(rd):
                cov[rd.reference_start] += 1
                cov[rd.reference_end] -= 1
    
def qv(fname):
    infos = []
    for line in open(fname):
        its = line.split()
        depth = int(its[3])
        match = its[4].count('.') + its[4].count(',')

        infos.append([depth, match, len(its[4])-depth])

    infos.sort(key=lambda x: x[0])
    vinfos = infos[len(infos) // 10 : 9*len(infos) // 10]
    
    cor = 0
    rd_err = [0, 0]
    for d, m, e in vinfos:
        if m*2 > d:
            cor += 1
            rd_err[0] += d
            rd_err[1] += d - m + e
        
    print(cor / len(vinfos))
    print(rd_err[1] / rd_err[0])
            
            





def main(argv):
    parser = argparse.ArgumentParser(description="")
    
    parser.add_argument("command", type=str)
    parser.add_argument("bam", type=str)
    parser.add_argument("--test", type=str)
    try:
    
        args = parser.parse_args(argv)
        if args.command == "n50":
            print(n50(args.bam))
        elif args.command == "coverage":
            coverage(args.bam)
        elif args.command == "qv":
            qv(args.test)

    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)



if __name__ == "__main__":
    main(sys.argv[1:])