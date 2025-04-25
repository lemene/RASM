#!/usr/bin/env python3

import os, sys
import argparse
import traceback
import yaml

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
    
def main(argv):
    parser = argparse.ArgumentParser(description="detect misassembly")
    parser.add_argument("command", type=str)
    parser.add_argument("--config", type=str, default="")
    parser.add_argument("--bam", type=str)
    parser.add_argument("--asm", type=str)
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--work-dir", default=".", help="work directory to output results")
    parser.add_argument("--min-contig", default=20000, type=int)  
    parser.add_argument("--dump", type=str, help="file for storing intermediate state")

    def make_summary(args):
        smry = summary.Summary(args.bam, args.asm)
        if args.dump:
            if not os.path.exists(args.dump) or utils.is_file_newer([args.bam, args.asm], args.dump):
                smry.scan(args.threads, args.min_contig)
                smry.save(args.dump)
            else:
                smry.load(args.dump) 

        else:
            smry.scan(args.threads, args.min_contig)
        return smry

    try:
        args = parser.parse_args(argv)

        if args.command == "detect":
            smry = make_summary(args)
            if args.config == "":
                args.config = os.path.join(utils.prj_dir, "..", "config", "default.yaml")
            
            config = yaml.safe_load(open(args.config))
            for k, v in vars(args).items():     # TODO Need to remove the unwanted
                config[k] = v
            detector = dt.Detector(config, args.work_dir, args.bam, args.asm)
            detector.summary = smry
            detector.detect(args.threads, args.min_contig)
        elif args.command == "n50":
            print(n50(args.bam))
        elif args.command == "coverage":
            coverage(args.bam)
        elif args.command == "summary":
            smry = make_summary(args)
            smry.stat(args.threads)
        elif args.command == "test":
            pass

    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)


if __name__ == "__main__":
    main(sys.argv[1:])