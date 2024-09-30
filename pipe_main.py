#!/usr/bin/env python3
import argparse
import os, sys
import traceback
import pysam
import info_stats
import detector as dt

import utils

def main():
    parser = argparse.ArgumentParser(description="detect misassembly")
    #group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-t", "--threads", dest="threads", type=int, default=1)
    parser.add_argument("--work-dir", dest="work_dir", default="output", help="work directory to output results")

    parser.add_argument("--ref", dest="ref", required=True)
    parser.add_argument("--fastq", dest="fastq", help="fastq reads file")
    parser.add_argument("--bam", dest="bam", required=True)
    parser.add_argument("--log", default="lrmd.log")
    
    parser.add_argument("--min-contig", dest="min_contig", help="skip SV consensus process contig shorter than this, keep with raw", default=20000, type=int)   # 200000 调 1,000,000     5,000,000
    
    parser.add_argument("--config", type=str)

    try:
        
        args = parser.parse_args()
        
        detector = dt.Detector(args.config, args.work_dir)

        utils.enable_logging(os.path.join(args.work_dir, "lrmd.log"))
        utils.logger.info("Start my pipe to assembly")

        args.work_dir = os.path.abspath(args.work_dir)

        utils.make_dir(args.work_dir)
        utils.make_dir(os.path.join(args.work_dir, "step2_candidate_regions"))

        args.bam = os.path.abspath(args.bam)
        args.ref = os.path.abspath(args.ref)

        bam = pysam.AlignmentFile(args.bam, "rb")
        ctgs = [(ref, rlen) for ref, rlen in zip(bam.references, bam.lengths) if rlen >= args.min_contig]

        detector.run_find_pipe(args.ref, args.bam, ctgs, os.path.join(args.work_dir, "step2_candidate_regions"), args.threads)

    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)

if __name__ == "__main__":
    main()
