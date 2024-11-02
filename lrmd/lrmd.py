#!/usr/bin/env python3

import os, sys
import argparse
import traceback
import pysam

import pysam
import detector as dt
import numpy as np
import multiprocessing as mp
import summary

import utils

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
    



def test(fname):
    bam = pysam.AlignmentFile(fname, "rb")
    #  CP132235.1:1226600-1227600
    for read in bam.fetch("CP132235.1", 1226600, 1227600):
        print(read.qname, read.mapping_quality, read.is_unmapped, read.is_secondary, read.is_supplementary)
    # return read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < mapq



def main(argv):
    parser = argparse.ArgumentParser(description="detect misassembly")
    parser.add_argument("command", type=str)
    parser.add_argument("--config", type=str)
    parser.add_argument("--bam", type=str)
    parser.add_argument("--asm", type=str)
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--work-dir", default=".", help="work directory to output results")
    parser.add_argument("--min-contig", default=20000, type=int)  
    parser.add_argument("--dump", type=str, help="file for storing intermediate state" )

    try:
        args = parser.parse_args(argv)

        if args.command == "detect":
            detector = dt.Detector(args.config, args.work_dir, args.bam, args.asm)
            detector.detect(args.threads, args.min_contig)
        elif args.command == "n50":
            print(n50(args.bam))
        elif args.command == "coverage":
            coverage(args.bam)
        elif args.command == "summary":
            smry = summary.Summary(args.bam, args.asm)
            smry.scan(args.threads, args.min_contig)
            smry.save("example.pkl")
            #smry.load("example.pkl") 
            #smry.stat(args.threads)
            #smry.stat_read()

        elif args.command == "test":
            import clip, pileup, depth
            smry = summary.Summary(args.bam, args.asm)
            smry.load("example.pkl") 
            
            bam = pysam.AlignmentFile(args.bam, "rb")
            ctgs = [(ref, rlen) for ref, rlen in zip(bam.references, bam.lengths) if rlen > args.min_contig]
            #info = clip.ClipInfo()
            #info.build(ctgs, args.bam, args.threads, 500)
            #info.get_mis_candidates(400, 200, 5)
            
            # info = pileup.PileupInfo()
            # info.build(ctgs, smry, args.threads, )
            # info.get_mis_candidates(400, 200, 0.2)

            info = depth.DepthInfo()
            info.build(ctgs, smry, args.threads )
            info.get_mis_candidates(400, 200, 0.15, 1.75)

    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)


if __name__ == "__main__":
    main(sys.argv[1:])