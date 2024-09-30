#!/usr/bin/env python

import pandas as pd
import numpy as np
import math
import sys,os
import pysam
import collections
import re
from multiprocessing import Pool
import shutil
import subprocess
import time

                                                            
            
def cluster_reg(ls, dis):
    if not ls: return []
    new_ls = []
    ctg = ls[0][0]
    start, end = -1, -1
    for reg in ls:
        if end > 0:
            if reg[1] - end > dis:
                new_ls.append([ctg, start, end])
                start, end = reg[1], reg[2]
            else:
                end = reg[2]
        else:
            start, end = reg[1], reg[2]
    if start > 0:
        new_ls.append([ctg, start, end])
    return new_ls

def write_pileup_dic(window_pileup_dict, f_out):
    data=pd.DataFrame(window_pileup_dict)
    data.to_csv(f_out, sep="\t", index=False)

def parse_reg_pileup(ctg, reg_start, reg_end, ref, bam, out_dir, params, data_type):
    ''' cal pileup data of a window'''
    t1 = time.time()
    
    
    win_size, step_size, cluster_dis = params["win_size"], params["step_size"], params["cluster_dis"]
    if data_type == "ont":
        min_correct_portion, max_differ_portion, max_disagree_portion = params['ont']["min_correct_portion"], params['ont']["max_differ_portion"], params['ont']["max_disagree_portion"]
    elif data_type == 'hifi':
        min_correct_portion, max_differ_portion, max_disagree_portion = params['hifi']["min_correct_portion"], params['hifi']["max_differ_portion"], params['hifi']["max_disagree_portion"]

    print("{}:{}-{}".format(ctg, reg_start, reg_end))
    region = ctg + ":" + str(reg_start) + "-" + str(reg_end - 1)
    # pileup_file = os.path.join(out_dir, region + "pileup.txt")
    # cmd = ["samtools mpileup", "-B", "-q", "20", "-aa", "-d 100","-r", region, "-f", ref, bam, ">", pileup_file]
    # subprocess.check_call(" ".join(cmd), shell=True)
    pileup_stream = pysam.mpileup("-B", "-q", "20", "-aa", "-d 100","-r", region, "-f", ref, bam)    # 计算pileup信息
    # print("cost:{}".format(time.time()- t1))
    # pileup_dict={"contig":[],"correct":[],"ambiguous":[],"insert":[],"deletion":[],"disagree":[],"depth":[]}
    pileup_dict={"contig":[],"correct":[],"ambiguous":[],"insert":[],"deletion":[],
                 "disagree":[],"depth":[],"differ":[]}
    for line in pileup_stream.split("\n"):
        record = line.split()
        # print(record)
        if not record:
            continue
        ## 
        pileup_dict['contig'].append(record[0])
        match_detail=record[4]
        pileup_dict['correct'].append(match_detail.count('.')+match_detail.count(','))
        # pileup_dict['ambiguous'].append(match_detail.count('*'))
        # pileup_dict['insert'].append(match_detail.count("+"))   # 只能记录插入的次数，反映不了插入大小
        # pileup_dict['deletion'].append(match_detail.count("-")) # 记录删除的次数
        pileup_dict["depth"].append(int(record[3]))     # 反映的是pileup计算的depth
        st = ''.join(re.split('[\+|\-][0-9]+[ATCGatcg]+', match_detail))
        disagree_numd = st.count('a')+st.count('A')+st.count('t')+st.count('T')+st.count('c')+st.count('C')+st.count('g')+st.count('G')  # 单碱基的不一致，snp
        pileup_dict["disagree"].append(disagree_numd)
        differ_numd = match_detail.count('a')+match_detail.count('A')+match_detail.count('t')+match_detail.count('T')+match_detail.count('c')+match_detail.count('C')+match_detail.count('g')+match_detail.count('G')  # 包括插入删除的大小信息
        pileup_dict["differ"].append(differ_numd)
    ## 按window计算信息
    ls = []
    # window_pileup_dict={"contig":[],"start_pos":[],"correct_portion":[],"ambiguous_portion":[],
    #                     "disagree_portion":[],"deletion_portion":[],"insert_portion":[],
    #                     "normalized_coverage":[],"normalized_deviation":[],"mean_coverage":[],
    #                     "differ_portion":[]}
    window_pileup_dict={"contig":[],"start_pos":[],"end_pos":[],"correct_portion":[],
                        "disagree_portion":[],"differ_portion":[]}
    for i in range(0, len(pileup_dict["correct"]), step_size):
        start = i
        end = i + win_size if i + win_size <= len(pileup_dict["correct"]) else len(pileup_dict["correct"])
        ## 
        win_start = start + reg_start
        win_end = win_start + win_size if win_start + win_size < reg_end else reg_end
        total = np.sum(pileup_dict["depth"][start:end])
        window_pileup_dict["contig"].append(ctg)
        window_pileup_dict["start_pos"].append(win_start)
        window_pileup_dict["end_pos"].append(win_end)
        if total == 0:
            window_pileup_dict["correct_portion"].append(np.nan)
            window_pileup_dict["differ_portion"].append(np.nan)
            window_pileup_dict["disagree_portion"].append(np.nan)
            ls.append([ctg, start + reg_start, end + reg_start])
            continue
        window_pileup_dict["correct_portion"].append(np.sum(pileup_dict['correct'][start:end])/total)
        # window_pileup_dict["ambiguous_portion"].append(np.sum(pileup_dict["ambiguous"][start:end])/total)        
        # window_pileup_dict["insert_portion"].append(np.sum(pileup_dict['insert'][start:end])/total)    
        # window_pileup_dict["deletion_portion"].append(np.sum(pileup_dict['deletion'][start:end])/total) 
        window_pileup_dict["disagree_portion"].append(np.sum(pileup_dict["disagree"][start:end])/total)
        window_pileup_dict["differ_portion"].append(np.sum(pileup_dict["differ"][start:end])/total) 
        # print("{}-{}: correct:{},differ:{}".format(start + reg_start, end + reg_start, window_pileup_dict["correct_portion"][-1], window_pileup_dict["differ_portion"][-1]))
        ## 
        if total == 0 \
            or window_pileup_dict["correct_portion"][-1] < min_correct_portion \
            or window_pileup_dict["differ_portion"][-1] > max_differ_portion \
            or window_pileup_dict["disagree_portion"][-1] > max_disagree_portion:
            ls.append([ctg, start + reg_start, end + reg_start])
    ls = cluster_reg(ls, cluster_dis)
    ## 
    bed_out = out_dir + "/" + ctg + ":" + str(reg_start) + "-" + str(reg_end) + "_pileup.bed"
    with open(bed_out, "w") as f:
        for reg in ls:
            f.write("{}\t{}\t{}\t{}\n".format(reg[0], reg[1], reg[2], str(reg[2]-reg[1])+"bp-pileup_reg"))
    pileup_out = out_dir + "/" + ctg + ":" + str(reg_start) + "-" + str(reg_end) + "_pileup_feature.txt"
    write_pileup_dic(window_pileup_dict, pileup_out)
    # return ls, window_pileup_dict


def parse_pileup_parallel(ctg_ls, ref, bam_in, threads, params, out_dir, data_type):
    # 
    print("Start pileup parse!!!")
    print("pileup_params: {}".format(params))
    reg_size = params["reg_size"]
    reg_out_dir = os.path.join(out_dir, "parts")
    if os.path.exists(reg_out_dir):shutil.rmtree(reg_out_dir)     # 目录存在要先删除，注意目录要给对，最好不是已有目录
    if not os.path.exists(reg_out_dir):os.makedirs(reg_out_dir)
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in + ".bai")
    # results = []
    
    pool = Pool(processes=threads)
    for ctg in ctg_ls:
        ctg_len = bam_reader.get_reference_length(ctg)
        print(ctg, "chr_len:", ctg_len)
        ## split ctg to reg
        for i in range(0, ctg_len, reg_size):
            reg_start = i
            reg_end = i + reg_size if i + reg_size <= ctg_len else ctg_len
            # task_ls.append([ctg, reg_start, reg_end])
            # print("{}:{}-{}".format(ctg, reg_start, reg_end))
            pool.apply_async(parse_reg_pileup, args=(ctg, reg_start, reg_end, ref, bam_in, reg_out_dir, params, data_type))
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    ## 合并所有数据
    bed_merge_cmd = ["cat", reg_out_dir+"/*.bed", "|", "sort -k 1,1 -k 2n,2", ">", out_dir+"/"+"candidate_pileup.bed"]
    subprocess.check_call(" ".join(bed_merge_cmd), shell=True)
    print("Run pileup parse Done!!!")


if __name__=='__main__':
    import yaml
    ctg_ls = sys.argv[1]
    ref = sys.argv[2]
    bam_in = sys.argv[3]
    threads = sys.argv[4]
    config_file = sys.argv[5]
    out_dir = sys.argv[6]

    ## process
    ctg_ls = ctg_ls.split(",")
    threads = int(threads)
    ref = os.path.abspath(ref)
    config_file = os.path.abspath(config_file)
    bam_in = os.path.abspath(bam_in)
    out_dir = os.path.abspath(out_dir)
    with open(config_file, "r") as f:
        config = yaml.safe_load(f.read())
        params = config["pileup_params"]
        # print("pileup_params: ", params)
    parse_pileup_parallel(ctg_ls, ref, bam_in, threads, params, out_dir)


