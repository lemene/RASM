#!/usr/bin/env python

import pandas as pd
import numpy as np
import math
import sys,os,argparse,warnings
import pysam
import collections
import re
from multiprocessing import Pool
import shutil
import subprocess
import time
def contig_pool(samfile):
    contig_len={}
    for (ref,lens) in zip(samfile.references,samfile.lengths):
        contig_len[ref]=lens
    return contig_len
         
def pileup_window_cal(pileup_dict):
    window_dict={"contig":[],"start_pos":[],"correct_portion":[],"ambiguous_portion":[],"disagree_portion":[],
    "deletion_portion":[],"insert_portion":[],"coverage":[],"deviation":[]}
    win_size = 1000
    # for i in range(300,len(pileup_dict['correct']),100):
    for i in range(0,len(pileup_dict['correct']),win_size):
        start=i
        end=i+1000
        # end = i + 300
        total=np.sum(pileup_dict['depth'][start:end])
        window_dict["contig"].append(pileup_dict["contig"][0])
        window_dict["start_pos"].append(start)
        window_dict["correct_portion"].append(np.sum(pileup_dict['correct'][start:end])/total)
        window_dict["ambiguous_portion"].append(np.sum(pileup_dict["ambiguous"][start:end])/total)        
        window_dict["insert_portion"].append(np.sum(pileup_dict['insert'][start:end])/total)    
        window_dict["deletion_portion"].append(np.sum(pileup_dict['deletion'][start:end])/total)   
        window_dict["disagree_portion"].append(np.sum(pileup_dict['disagree'][start:end])/total) 
        window_dict["coverage"].append(np.mean(pileup_dict["depth"][start:end]))
        window_dict["deviation"].append(np.sqrt(np.var(pileup_dict["depth"][start:end]))/np.mean(pileup_dict["depth"][start:end]))        
        if len(pileup_dict['correct']) - (i+100) <= 300:
            break
    return window_dict                                                        
            
def pileupfile_parse(args):
    """
    process pileup file
    """    
    samfile=pysam.AlignmentFile(args.bam,"rb") 
    contig_len=contig_pool(samfile)
    prev_contig=None    
    pileup_dict={"contig":[],"correct":[],"ambiguous":[],"insert":[],
        "deletion":[],"disagree":[],"depth":[]} 
    window_pileup_dict={"contig":[],"start_pos":[],"correct_portion":[],"ambiguous_portion":[],"disagree_portion":[],
    "deletion_portion":[],"insert_portion":[],"normalized_coverage":[],"normalized_deviation":[],"mean_coverage":[]}  
    for line in open(args.pileup,"r"):
        record = line.strip().split('\t')
        if contig_len[record[0]] < args.mlen:
            continue
        if prev_contig is None:
            prev_contig=record[0]
        if record[0] !=prev_contig:
            window_data=pileup_window_cal(pileup_dict)
            mean_cov=np.mean(window_data["coverage"])
            window_pileup_dict["contig"].extend(window_data["contig"]) 
            window_pileup_dict["start_pos"].extend(window_data["start_pos"]) 
            window_pileup_dict["correct_portion"].extend(window_data["correct_portion"]) 
            window_pileup_dict["ambiguous_portion"].extend(window_data["ambiguous_portion"]) 
            window_pileup_dict["disagree_portion"].extend(window_data["disagree_portion"]) 
            window_pileup_dict["deletion_portion"].extend(window_data["deletion_portion"])   
            window_pileup_dict["insert_portion"].extend(window_data["insert_portion"])     
            window_pileup_dict["normalized_coverage"].extend(window_data["coverage"]/mean_cov)
            window_pileup_dict["normalized_deviation"].extend(window_data["deviation"])
            window_pileup_dict["mean_coverage"].extend([mean_cov]*len(window_data["start_pos"]))                                           
            pileup_dict={"contig":[],"correct":[],"ambiguous":[],"insert":[],"deletion":[],"disagree":[],"depth":[]}
            prev_contig = record[0]
        pileup_dict['contig'].append(record[0])
        match_detail=record[4]        
        pileup_dict['correct'].append(match_detail.count('.')+match_detail.count(','))
        # pileup_dict['ambiguous'].append(match_detail.count('*'))
        # pileup_dict['insert'].append(match_detail.count("+"))
        # pileup_dict['deletion'].append(match_detail.count("-"))
        pileup_dict['ambiguous'].append(match_detail.count('*'))
        pileup_dict['insert'].append(match_detail.count("+"))   # 只能记录插入未知的数目，反映不了插入具体情况
        pileup_dict['deletion'].append(match_detail.count("-"))
        
        pileup_dict['depth'].append(int(record[3]))
        st = ''.join(re.split('[\+|\-][0-9]+[ATCGatcg]+',match_detail))
        numd = st.count('a')+st.count('A')+st.count('t')+st.count('T')+st.count('c')+st.count('C')+st.count('g')+st.count('G')
        pileup_dict['disagree'].append(numd)
        # print(pileup_dict)
    ## 
    window_data=pileup_window_cal(pileup_dict)
    mean_cov=np.mean(window_data["coverage"])   # cal mean cov of the contig
    window_pileup_dict["contig"].extend(window_data["contig"])
    window_pileup_dict["start_pos"].extend(window_data["start_pos"])
    window_pileup_dict["correct_portion"].extend(window_data["correct_portion"])
    window_pileup_dict["ambiguous_portion"].extend(window_data["ambiguous_portion"])
    window_pileup_dict["disagree_portion"].extend(window_data["disagree_portion"])
    window_pileup_dict["deletion_portion"].extend(window_data["deletion_portion"])
    window_pileup_dict["insert_portion"].extend(window_data["insert_portion"])     
    window_pileup_dict["normalized_coverage"].extend(window_data["coverage"]/mean_cov)
    window_pileup_dict["normalized_deviation"].extend(window_data["deviation"])
    window_pileup_dict["mean_coverage"].extend([mean_cov]*len(window_data["start_pos"]))

    if not os.path.exists(os.path.join(args.output, "temp/pileup")):
        os.makedirs(os.path.join(args.output, "temp/pileup"))
    data=pd.DataFrame(window_pileup_dict)
    data.to_csv(os.path.join(args.output, "temp/pileup/pileup_feature.txt"),sep="\t")
    ## 
    return data

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

def call_back(res):
    # return
    print(res)
def error_call_back(error_code):
    # return
    print("error: ", error_code)

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

def test():
    ref = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe2/corrected_ref/reference.fasta"
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe2/step1_mapping/aln.sorted.bam"
    ctg_ls = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in + ".bai").references
    threads = 10
    params = {"win_size":1000, "step_size":500, "min_correct_portion":0.9, "max_differ_portion":0.1, "max_disagree_portion":0.02, "cluster_dis":1000, "reg_size":1000000}
    out_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe2/step2_candidate_regions/pileup"
    parse_pileup_parallel(ctg_ls, ref, bam_in, threads, params, out_dir)
def test1():
    ref = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe2/corrected_ref/reference.fasta"
    bam = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe2/step1_mapping/aln.sorted.bam"
    ctg = "chrX"
    start = 0
    end = 725652
    params = {"win_size":1000, "step_size":500, "min_correct_portion":0.9, "max_differ_portion":0.1, "max_disagree_portion":0.02, "cluster_dis":1000, "reg_size":1000000}
    out_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/yeast_ont/my_pipe2/step2_candidate_regions/pileup/parts"
    parse_reg_pileup(ctg, start, end, ref, bam, out_dir, params)

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
    # test()
    # test1()
    # import time
    # t1 = time.time()
    # ctg = "NC_000023.11"
    # # start = 115000000
    # # end = 116000000
    # start = 135000000
    # end = 136000000
    # # start = 135706000
    # # end = 135714000
    # # start = 0
    # # end = 20000
    # ref = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/corrected_ref/reference.fasta"
    # bam = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/step1_mapping/aln.sorted.bam"
    # params = {"win_size":1000, "step_size":500, "min_correct_portion":0.9, "max_differ_portion":0.1, "max_disagree_portion":0.02, "cluster_dis":1000, "reg_size":1000000}
    # out_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/pileup/parts"
    # # ls, dic = parse_reg_pileup(ctg, start, end, ref, bam, 200, 100, 0.9, 0.1, 1000)
    # parse_reg_pileup(ctg, start, end, ref, bam, out_dir, params)
    # print(time.time() - t1)

    # ## 
    # # ctg_ls = ["NC_000023.11"]
    # ref = ref
    # bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/step1_mapping/aln.sorted.bam"
    # ctg_ls = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in + ".bai").references
    # threads = 40
    # params = params
    # out_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/pileup/all"
    # # parse_pileup_parallel(ctg_ls, ref, bam_in, threads, params, out_dir)
    # print(time.time() - t1)

    # ## 
    # pass

