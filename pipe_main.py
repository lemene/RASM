#!/usr/bin/env python3
import argparse
import os, sys
import traceback
import pysam
import yaml
from create_consensus_by_bed import fasta_parser, scaffold_by_rec, Utils
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
    
    ## 
    parser.add_argument("--min-contig", dest="min_contig", help="skip SV consensus process contig shorter than this, keep with raw", default=20000, type=int)   # 200000 调 1,000,000     5,000,000
    
    parser.add_argument("--config", type=str)

    try:
        
        args = parser.parse_args()
        
        detector = dt.Detector()

        utils.logger.info("Start my pipe to assembly")
        ### 绝对路径配置
        args.work_dir = os.path.abspath(args.work_dir)

        ### 生成文件路径配置
        work_dir = args.work_dir
        step2_dir = work_dir + "/" + "step2_candidate_regions"
        log_dir = work_dir + "/" + "logs"
        config_dir = work_dir + "/" + "config"
        for Dir in [work_dir, step2_dir, log_dir, config_dir]:
            Utils.make_dir(Dir)
        assembly_log = os.path.join(log_dir, "assembly.log")
        utils._enable_logging(assembly_log, debug=False, overwrite=True)

        args.bam = os.path.abspath(args.bam)
        args.ref = os.path.abspath(args.ref)

        ## 一些参数的获取
        config = yaml.safe_load(open(args.config))
        utils.logger.info("Config参数设置: {}".format(config))

        bam = pysam.AlignmentFile(args.bam, "rb")
        ctgs = [(ref, rlen) for ref, rlen in zip(bam.references, bam.lengths) if rlen >= args.min_contig]

        detector.run_find_pipe(args.ref, args.bam, ctgs, step2_dir, args.threads, config["step2"])

        info_stats.get_qv(step2_dir, ctgs, args.threads)

    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)

if __name__ == "__main__":
    main()
