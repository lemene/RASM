#! /usr/bin/env python3

import sys,os
import multiprocessing
import argparse
import traceback


def build_bam(ref, rds, bam, threads, params):
    cmd = "minimap2 %s %s -t %d %s -a | samtools sort -@ %d > %s" % (ref, rds, threads, params, threads, bam) 
    print(cmd)
    r = os.system(cmd)
    assert r == 0
    cmd = "samtools index -@ %d %s" % bam
    r = os.system(cmd)
    assert r == 0

def main(argv):
    parser = argparse.ArgumentParser("map reads to reference and generate bamfile")
    parser.add_argument("reference", type=str)
    parser.add_argument("reads", type=str)
    parser.add_argument("bam", type=str)
    parser.add_argument("threads", type=int)
    parser.add_argument("params", type=str)
    try:
        args = parser.parse_args(argv)
        build_bam(args.reference, args.reads, args.bam, args.threads, args.params)

    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)


if __name__ == '__main__':
    main(sys.argv[1:])

