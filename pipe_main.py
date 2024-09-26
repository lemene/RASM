#!/usr/bin/env python3
import argparse
import os
import sys
import pysam
# import signal
# import multiprocessing as mp
import time
import subprocess
import shutil
import logging
import yaml
from find_candidate_regions.find_candidate_regions import run_find_candidate2_parallel
from create_consensus_by_bed.get_fasta_consensus2 import run_SVconsensus_parallel
from create_consensus_by_bed import fasta_parser, scaffold_by_rec, Utils
from find_candidate_regions.find_mis_pipe import run_find_pipe
import info_stats
logger = logging.getLogger()
def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


def run_minimap2(target_fa, query_fa, data_type, threads, minimap_out, params_ls):     
    ''' 运行minimap2 提供params '''
    t1 = time.time()
    params = " ".join([str(param) for param in params_ls])
    minimap_cmd = ["/usr/bin/time -v minimap2", "-ax", "map-"+data_type, params, target_fa, query_fa, "-t", str(threads), ">", minimap_out]
    logger.info("Running: %s", " ".join(minimap_cmd))
    subprocess.check_call(" ".join(minimap_cmd), shell=True)
    logger.info("Run minimap2 finished, cost {}s".format(time.time() - t1))

def run_samtools(sam_in, samtools_out, threads):
    t1 = time.time()
    samtools_cmd = ["/usr/bin/time -v samtools", "sort", "-O", "BAM", "-@", str(threads), sam_in, "-o", samtools_out, "&&", "samtools", "index", "-@", str(threads), samtools_out]
    logger.info("Running: %s", " ".join(samtools_cmd))
    subprocess.check_call(" ".join(samtools_cmd), shell=True)
    logger.info("Run samtools finished, cost {}s".format(time.time() - t1))


def main():
    parser = argparse.ArgumentParser(description="reference guided asssembly")
    #####
    ### ,MIN_CLIP_NUM, MIN_DEP, MIN_DEP_REG
    parser = argparse.ArgumentParser(description="get_candidate_regions")
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-t", "--threads", dest="threads", type=int, default=1)
    group.add_argument("--ctg-ls", dest="ctg_ls", help="ctg list to perform assembly, recommend to provide contig ls, only accept format like:chr1,chr2")   # contig list，影响结果
    group.add_argument("-a", dest="all_chrs", default=False, action='store_true', help="we will keep and process all chrs")
    ### 路径
    parser.add_argument("--work-dir", dest="work_dir", required=True, help="work directory to output results")
    # parser.add_argument("--bam", dest="bam_in", required=True, help="bam file")
    parser.add_argument("--ref", dest="ref", required=True)
    parser.add_argument("--fastq", dest="fastq", help="fastq reads file")
    ## 
    parser.add_argument("--data-type", dest="data_type", required=True, choices=["ont", "hifi", "pb"], help="fastq file type")
    ### 性能参数
    # parser.add_argument("--min_clip_num", dest="min_clip_num", default=5, help="min_clip_num of window to be selected") # important
    # parser.add_argument("--min_clip_len", dest="min_clip_len", default=500, help="min_clip_len of clip to be selected as clip")
    # parser.add_argument("--dp-win-size", dest="dp_win_size", default=100, help="dp win_size")
    # parser.add_argument("--min-dep", dest="min_dep", default=10, help="dep threshold to be select as low dep region")
    # parser.add_argument("--max-dep", dest="max_dep", default=48, help="max dep threshold")  # ->后面改为随由平均cov来计算 (1.7 * avg cov)
    # parser.add_argument("--min-dep-reg", dest="min_dep_reg", default=100, help="minimum length of low dep region")
    parser.add_argument("--min-contig", dest="min_contig", help="skip SV consensus process contig shorter than this, keep with raw", default=20000, type=int)   # 200000 调 1,000,000     5,000,000
    parser.add_argument("--fill-size", dest="Nfill_size", type=int, default=100, help="The N filling size for the assembly gap region.")
    # 区间聚类和延展参数
    # parser.add_argument("--cluster-params", dest="cluster_params")
    ### 选择性的性能参数，
    parser.add_argument("--ex-unmapped-denovo", dest="ex_unmapped_denovo", default=False, action="store_true", help="exclude unmapped denovo asm, may be unique sequence of sample")    # 默认保留，获得更加完整的组装结果
    '''Test'''
    # parser.add_argument("--cut-gap", dest="cut_gap", default=False, action="store_true", help="cut between big gap instead of N_fill")      
    parser.add_argument("--polish-filter", dest="polish_filter", default=False, action="store_true")
    ## 
    parser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true", help="overwrites all existing results")
    parser.add_argument("-R", dest="re_run", choices=["step1", "step2", "step3", "step4", "step5"], default=False, help="re_run step")
    ## 1、mis find mode
    # parser.add_argument("--min_span", dest="min_span", type=int, default=5)
    # parser.add_argument("--boundary", dest="boundary", default=3000, type=int)
    parser.add_argument("--bam", dest="bam")
    ## 2、correct mode。纠错模式？？？
    parser.add_argument("--correct", dest="correct", default=False, action="store_true")
    parser.add_argument("--min_corr_len", dest="min_corr_len", type=int, default=500000)
    parser.add_argument("--skip_denovo", dest="skip_denovo", default=False, action="store_true")
    parser.add_argument("--cut", dest="cut", default=False, action="store_true")
    # parser.add_argument("")
    ## 一些参数的配置文件
    parser.add_argument("--config", default='/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Configs/Config.yaml')
    ## 
    parser.add_argument("--out_denovo", dest="out_denovo", default=False, help="Provide path to the denovo fa, for testing")

    ##
    args = parser.parse_args()
    

    logger.info("Start my pipe to assembly")
    ### 绝对路径配置
    args.work_dir = os.path.abspath(args.work_dir)
    ### 生成文件路径配置
    work_dir = args.work_dir
    ref_dir = work_dir + "/" + "corrected_ref"
    step1_dir = work_dir + "/" + "step1_mapping"
    step2_dir = work_dir + "/" + "step2_candidate_regions"
    log_dir = work_dir + "/" + "logs"
    config_dir = work_dir + "/" + "config"
    dp_file_dir = work_dir + "/" + "depths"
    for Dir in [work_dir, ref_dir, dp_file_dir, step1_dir, step2_dir, log_dir, config_dir]:
        Utils.make_dir(Dir)
    ## log和文件检查
    assembly_log = os.path.join(log_dir, "assembly.log")
    _enable_logging(assembly_log, debug=False, overwrite=True)

    args.bam = os.path.abspath(args.bam)
    args.ref = os.path.abspath(args.ref)
    ## 一些参数的获取
    with open(args.config, "r") as f: # config参数获取
        config = yaml.safe_load(f.read())     # 获取部分参数
        logger.info("Config参数设置: {}".format(config))
    config_save = os.path.join(config_dir, "WorkConfig.yaml")
    shutil.copy(args.config, config_save)
    
    ## 参考序列的处理
    corrected_ref = ref_dir + "/" + "reference.fasta"

    print("Find mis mode, Skip reference process")
    if not os.path.isfile(args.ref + ".fai"):
        pysam.faidx(args.ref)
    # shutil.copy(args.ref, corrected_ref)
    ctg_ls = []
    with open(args.ref + ".fai", "r") as f:
        for line in f:
            chr = line.strip().split("\t")[0]
            ctg_ls.append(chr)

    ## 
    overwrite = args.overwrite
    re_run = args.re_run
    if re_run: print("Will re run from {}".format(re_run))

######################################################### pipeline #########################################################
    
    ## step1 mapping 
    t1 = time.time()
    print("Find mis mode, Skip mapping")
    samtools_out1 = args.bam
    overwrite = True
    
    ## step2 Find candidate regions
    t2 = time.time()
    print("STEP1 cost: ", t2 - t1)

    keep_ls = []
    bam_reader = pysam.AlignmentFile(samtools_out1, "rb", index_filename=samtools_out1+".bai")
    print("--------correct/find mis mode-------")
    min_contig = args.min_contig
    new_ctg_ls = []
    for chr in ctg_ls:
        if bam_reader.get_reference_length(chr) > min_contig:
            new_ctg_ls.append(chr)
        else:
            keep_ls.append(chr)
    print("We will corr/find mis on:{}".format(new_ctg_ls))
    ctg_ls = new_ctg_ls

    # 1、mis find1
    run_find_pipe(args.ref, samtools_out1, ctg_ls, step2_dir, args.threads, config["step2"])
    
    # 2、some stats
    print("Start stats")
    ctg_ls = info_stats.simple_stats(args.ref, min_contig, step2_dir)
    info_stats.get_qv(args.ref, args.bam, step2_dir, ctg_ls, args.threads)


if __name__ == "__main__":
    main()
