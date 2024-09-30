import yaml
import os
import math

import multiprocessing as mp

import utils

from find_mis_pipe import *
import os
import numpy as np
import gzip
from multiprocessing import Pool
import subprocess
import pysam
import time
from collections import namedtuple
chr_info = namedtuple('chr_info', ["chr_id", "chr_len"])
Region = namedtuple('Region', ["chr_id", "start", "end"])

######################## Depth parse ##############################



def run_mosdepth2(ctg, bam_in, out_dir, DP_WIN_SIZE, threads, min_MQ):
    # min_MQ = 20
    # DP_WIN_SIZE = 100
    prefix = out_dir+"/"+ctg
    ## -n -x 用于缩短时间
    cmd = ["mosdepth", "-Q", str(min_MQ), "-b", str(DP_WIN_SIZE), "-t", str(threads), "-c", ctg, "-n", prefix, bam_in]    # mosdepth -b 100 test/NC_19 ../step1_mapping/aln.sorted.BAM      # 指定winsize
    print("Running: %s", " ".join(cmd))
    subprocess.check_call(" ".join(cmd), shell=True)
    print("Run {} done".format(ctg))

def read_mosdepth_dp_file(dp_file):
    '''best for ctg dp file'''
    dp_ls = []
    with gzip.open(dp_file, "rt") as f:
        lines = f.readlines()
        for line in lines:
            fields = line.strip().split("\t")
            dp = float(fields[3])
            dp_ls.append(dp)
    return dp_ls


def get_dp(ctg, bam_in, ctg_out_dir, DP_WIN_SIZE, min_MQ):
    '''per ctg'''
    if not os.path.isdir(ctg_out_dir):os.makedirs(ctg_out_dir)
    dp_out_file = ctg_out_dir+"/"+ctg + ".regions.bed.gz"
    ## 
    run_mosdepth2(ctg, bam_in, ctg_out_dir, DP_WIN_SIZE, 4, min_MQ)
    dp_ls = read_mosdepth_dp_file(dp_out_file)
    return {ctg: dp_ls}     # 返回一个字典存储

def get_dp_info_parallel(ctgs, bam_in, threads, out_dir, DP_WIN_SIZE, Block_size, min_MQ):
    print("----------------get_dp_info_parallel----------------")
    dp_file_dir = os.path.join(out_dir, "depths")
    dp_info_dir = os.path.join(out_dir, "dp_info")
    if not os.path.isdir(dp_file_dir):os.makedirs(dp_file_dir)
    if not os.path.isdir(dp_info_dir):os.makedirs(dp_info_dir)

    dp_dic = {}
    pool = Pool(processes=threads)
    results = [pool.apply_async(get_dp, args=(ctg[0], bam_in, dp_file_dir + "/" + ctg[0], DP_WIN_SIZE, min_MQ)) for ctg in ctgs]
    pool.close() 
    pool.join()

    for res in results:
        dp_dic.update(res.get())    # 更新所有的

    whole_dp_ls = []
    for ctg, ctg_dp_ls in dp_dic.items():
        whole_dp_ls.extend(ctg_dp_ls)
    whole_dp = np.median(whole_dp_ls)   # 全局的dp，使用中位数表示
    print("Whole dp: {} !!!".format(whole_dp))

    ## 
    dpinfo_dic = {}
    
    for ctg, ctg_len in ctgs:
        dpinfo_dic[ctg] = Depth_info(ctg, ctg_len, ctg_dp_ls, DP_WIN_SIZE, Block_size, whole_dp)
    
    return dpinfo_dic


def cluster_by_dis(reg_ls_in, dis): # 根据距离进行聚类
    if len(reg_ls_in) <= 1:
        return reg_ls_in
    reg_ls_in = sorted(reg_ls_in, key=lambda x: x[1])   # 聚类之前先排序
    reg_start = -1
    reg_end = -1
    chr_id = reg_ls_in[0].chr_id
    reg_ls_out = []
    for reg in reg_ls_in:
        if reg_start > -1:  # 非首次
            if reg.start - reg_end <= dis:
                reg_end = reg.end
                need_to_cluster.append(reg)
            else:   # new_reg
                reg_ls_out.append(Region(chr_id, reg_start, reg_end))
                reg_start = reg.start
                reg_end = reg.end
                need_to_cluster = [reg]
        else:
            reg_start = reg.start
            reg_end = reg.end
            need_to_cluster = [reg]
    if reg_start > -1:
        reg_ls_out.append(Region(chr_id, reg_start, reg_end))
    return reg_ls_out

def find_by_dp(dp_ls, dp_win_size, CHR_INFO:chr_info, MIN_DP, MAX_DP, bed_out):  # 传如dp列表，以及dp计算的window大小
    '''
    find depth reg by cov, 
    
    '''
    length = len(dp_ls)
    chr_len = int(CHR_INFO.chr_len)
    chr_id = CHR_INFO.chr_id

    ls = []
    # pre_dp = 0
    for i in range(length):
        if dp_win_size * (i + 1) < chr_len:
            if dp_ls[i] < MIN_DP or dp_ls[i] > MAX_DP:
                ls.append(Region(chr_id, dp_win_size * i, dp_win_size * (i + 1)))
        else:
            if dp_ls[i] < MIN_DP or dp_ls[i] > MAX_DP:
                ls.append(Region(chr_id, dp_win_size * i, chr_len))
    final_ls = cluster_by_dis(ls, 1000)
    with open(bed_out, "w") as fout:
        for reg in final_ls:
            fout.write("{}\t{}\t{}\t{}bp-cov_reg\n".format(reg.chr_id, reg.start, reg.end, reg.end-reg.start))
            # print("{}\t{}\t{}".format(reg.chr_id, reg.start, reg.end))
    return final_ls

class Detector:
    def __init__(self, cfg, wrkdir):
        self.cfg = yaml.safe_load(open(cfg))
        self.wrkdir = wrkdir

        utils.make_dir(wrkdir)

        utils.logger.info("Config参数设置: {}".format(self.cfg))

    def make_dp_info():
        pass


    def run_find_pipe(self, ref, bam, ctgs, out_dir, threads):
        config = self.cfg["step2"]

        dp_win_size = config["dp_params"]["dp_win_size"]
        block_size = config["dp_params"]["block_size"]
        min_MQ = config['dp_params']['min_MQ']
        dp_info_dic = get_dp_info_parallel(ctgs, bam, threads, out_dir, dp_win_size, block_size, min_MQ)

        ## ----------------1、find candidate----------------
        # dirs
        candidate_dir = os.path.join(out_dir, "candidate")
        if not os.path.isdir(candidate_dir):os.mkdir(candidate_dir)
        info_bed_dir = os.path.join(out_dir, "info")
        if not os.path.isdir(info_bed_dir):os.mkdir(info_bed_dir)
        filtered2_dir = os.path.join(out_dir, "filtered2")
        if not os.path.isdir(filtered2_dir):os.mkdir(filtered2_dir)
        pileup_dir = os.path.join(out_dir, "pileup")
        if not os.path.isdir(pileup_dir):os.mkdir(pileup_dir)

        # run find
        pool = Pool(processes=threads)
        for ctg in ctgs:
            pool.apply_async(find_candidate, args=(ctg[0], 0, ctg[1], ref, bam, out_dir, dp_info_dic[ctg[0]], config))
        
        pool.close() 
        pool.join()



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