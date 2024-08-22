import pysam
import sys
import os
import re
import random
import signal
import multiprocessing as mp
from multiprocessing import Pool
import argparse
from collections import namedtuple
import gzip
import subprocess
import pandas as pd
# from find_reg_by_depth import find_by_dp, find_by_dp2, get_dp_info_parallel, Depth_info      # vscode调
from find_candidate_regions.find_reg_by_depth import find_by_dp, find_by_dp2, get_dp_info_parallel, Depth_info    # shell调
from find_candidate_regions.find_from_pileup import parse_pileup_parallel
Region = namedtuple('Region', ["chr_id", "start", "end"])   
chr_info = namedtuple('chr_info', ["chr_id", "chr_len"])

def cal_sum(arry, l, r):
    # ans = 0
    # for i in range(l, r):
    #     ans += arry[i]
    return sum(arry[l:r])

def get_ctg_len(fai, ctg):
    ctg_len = 0
    with open(fai, "r") as f:
        for line in f:
            if line.split("\t")[0] == ctg:
                ctg_len = int(line.split("\t")[1])
    return ctg_len

def cal_depth(bam_in, ctg):  # cal depth of a contig, return list of depth, length == contig length
    bam_index = bam_in + ".bai"
    '''time samtools depth -a -J -Q $min_MQ -r $REG ${bam_file} > $out_prefix.depth    
    # del处depth不为0 -J会把deletion算到depth里面去(在deletion处depth不为0)'''
    min_MQ = 20
    min_len = 2000
    ctg_len = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index).get_reference_length(ctg)
    depth_ls = [0] * ctg_len
    # depth_stream = pysam.depth("-a", "-J", "-Q", str(min_MQ), "-r", ctg, bam_in)    # depth_stream接收stdout
    depth_stream = pysam.depth("-a", "-J", "-Q", str(min_MQ), "-l", str(min_len), "-g", "0x800", "-r", ctg, bam_in)
    for line in depth_stream.split("\n"):
        line_ls = line.split()
        if not line_ls:
            continue
        depth_ls[int(line_ls[1]) - 1] = int(line_ls[2])    # chr_id pos depth   1-based
    return depth_ls

def cluster_regions(reg_in):    # for clip regioins   return clustered regions in list type    input just region_ls of a contig
    if len(reg_in) <= 1:
        return reg_in
    ctg = reg_in[0].chr_id
    print("cluster regions of {} !!!".format(ctg))
    reg_start = -1
    reg_end = -1
    reg_out = []
    max_len = 1000   # 500 1000 1500
    need_to_cluster = []
    for reg in reg_in:
        # if reg[0] == reg_end:   # 区间有连接
        #     reg_end = reg[1]
        if reg_start > -1:  # 非首次
            if reg.start - reg_end <= max_len:
                reg_end = reg.end
                need_to_cluster.append(reg)
            else:   # new_reg
                # if reg_start > -1:   # 非首次
                #     reg_out.append(Region(ctg, reg_start, reg_end))
                    # reg_start = reg.start
                    # reg_end = reg.end
                # print("cluster from {} \n-> {}".format(need_to_cluster, Region(ctg, reg_start, reg_end)))
                reg_out.append(Region(ctg, reg_start, reg_end))
                reg_start = reg.start
                reg_end = reg.end
                need_to_cluster = [reg]
        else:
            reg_start = reg.start
            reg_end = reg.end
            need_to_cluster = [reg]
    if reg_start > -1:
        reg_out.append(Region(ctg, reg_start, reg_end))
    return reg_out

def cluster_regions2(reg_in):    # lowdep cluster return clustered regions in list type
    if len(reg_in) < 1:
        return []
    ctg = reg_in[0].chr_id
    reg_start = -1
    reg_end = -1
    reg_out = []
    max_len = 200
    for reg in reg_in:
        # if reg[0] == reg_end:   # 区间有连接
        #     reg_end = reg[1]
        if reg_start > -1:
            if reg.start - reg_end <= max_len:
                reg_end = reg.end
            else:   # new_reg
                # if reg_start > -1:   # 非首次
                #     reg_out.append(Region(ctg, reg_start, reg_end))
                    # reg_start = reg.start
                    # reg_end = reg.end
                reg_out.append(Region(ctg, reg_start, reg_end))
                reg_start = reg.start
                reg_end = reg.end
        else:
            reg_start = reg.start
            reg_end = reg.end
    if reg_start > -1:
        reg_out.append(Region(ctg, reg_start, reg_end))
    return reg_out

def find_from_clip(bam_in, ctg, bed_out, MIN_CLIP_NUM, MIN_CLIP_LEN):  # find candidate regions by clips
    # region_ls_merge = []    # 保存所有ctg的region
    bam_index = bam_in + ".bai"
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index)
    ctg_len = bam_reader.get_reference_length(ctg)  # 172126628
    reg_start = 0
    reg_end = ctg_len

    print("{}:{}-{}".format(ctg, reg_start, reg_end))
    clip_ls = [0] * ctg_len # 记录contig每个位置>min_clip_len的clip数目
    MIN_CLIP_LEN = 1000      # 500 300
    MIN_MAPPING_QUALITY = 0    # 10
    for read in bam_reader.fetch(ctg, start=reg_start, end=reg_end):
        if read.is_unmapped or read.is_secondary or read.mapping_quality < MIN_MAPPING_QUALITY:
            continue
        ref_pos = read.reference_start
        cigar = read.cigarstring
        tokens = re.findall("[\d]{0,}[A-Z]{1}", cigar)
        # left
        left = tokens[0]
        if left[-1] in "HS":
            if int(left[:-1]) > MIN_CLIP_LEN:
                # clip_ls.append([read.reference_start, int(left[:-1])]) # 记录下clip的位置，长度
                clip_ls[read.reference_start] += 1
        
        # right
        right = tokens[-1]
        if right[-1] in "HS":
            if int(right[:-1]) > MIN_CLIP_LEN:
                # clip_ls.append([read.reference_end, int(left[:-1])]) # 记录下clip的位置，长度
                # if read.reference_end > len(clip_ls)-1:
                #     print(read.query_name, "\t"+ctg+":"+str(read.reference_start)+"-"+str(read.reference_end))
                clip_ls[read.reference_end-1] += 1
    # clip_ls.sort(key=lambda clip:clip[0])
    print(ctg+":"+"FIND clip_regions")
    win_size = 1000  # 200 400
    stride = 500
    # MIN_CLIP_NUM = 5
    win_num = ctg_len // stride + 1 # 
    win_clip_num = [0] * win_num    # clip num of a window

    region_ls = []
    for i in range(win_num):    # i*stride, i*stride+win_size
        win_clip_num[i] = cal_sum(clip_ls, i*stride, i*stride+win_size) if i*stride+win_size < ctg_len else cal_sum(clip_ls, i*stride, ctg_len)
    
    for i in range(win_num):
        if win_clip_num[i] >= MIN_CLIP_NUM:
            if i*stride+win_size < ctg_len:
                # print("{}:{}-{}: {} clip".format(ctg, i*stride, i*stride+win_size, win_clip_num[i]))
                region_ls.append(Region(ctg, i*stride, i*stride+win_size))
            else:
                # print("{}:{}-{}: {} clip".format(ctg, i*stride, ctg_len, win_clip_num[i]))
                region_ls.append(Region(ctg, i*stride, ctg_len))
    region_ls = cluster_regions(region_ls)
    # region_ls_merge.extend(region_ls)

    # bed_out = os.path.join(out_dir, ctg, "clip.bed")
    fout = open(bed_out, "w")
    for region in region_ls:
        fout.write("{}\t{}\t{}\t{}\n".format(region.chr_id, region.start, region.end, str(region.end-region.start)+"bp-clip_reg"))
        # print("{}:{}-{} clip_reg".format(region.chr_id, region.start, region.end))
    print(ctg+":"+"FIND clip_regions Done!!!")
    fout.close

def bed_merge(fin_ls, fout):  # format: chr start   end   INFO
    f_out = open(fout, "w")
    out_ls = []
    for fin in fin_ls:
        with open(fin, "r") as f_in:
            for line in f_in:
                ctg, start, end, INFO = line.strip().split()[:4]
                # if INFO.endswith("low_dep"):
                #     f_out.write("{}\t{}\t{}\t{}\n".format(ctg, start, end, "low_dep"))
                # elif INFO.endswith("clip_reg"):
                #     f_out.write("{}\t{}\t{}\t{}\n".format(ctg, start, end, "clip_reg"))
                # else:
                #     f_out.write("{}\t{}\t{}\t{}\n".format(ctg, start, end, INFO))
                ##
                if INFO.endswith("low_dep"):
                    out_ls.append([ctg, start, end, "low_dep"])
                elif INFO.endswith("clip_reg"):
                    out_ls.append([ctg, start, end, "clip_reg"])
                elif INFO.endswith("cov_reg"):
                    out_ls.append([ctg, start, end, "clip_reg"])    # 懒得搞了
                elif INFO.endswith("pileup_reg"):
                    out_ls.append([ctg, start, end, "clip_reg"])    # 懒得搞了
                else:
                    out_ls.append([ctg, start, end, INFO])
    ## sort out_ls
    out_ls.sort(key=lambda res:(res[0], int(res[1])))
    for res in out_ls:
        ctg, start, end, info = res
        f_out.write("{}\t{}\t{}\t{}\n".format(ctg, start, end, info))
    f_out.close()
        

def get_dp_ls(dp_file):
    dp_ls = []
    with gzip.open(dp_file, "rt") as f:
        for line in f:
            fields = line.strip().split()
            dp_ls.append(float(fields[3]))
    return dp_ls

def run_find_candidate(bam_in, ctg_ls, out_dir, dp_file_dir, MIN_CLIP_NUM, MIN_CLIP_LEN, DP_WIN_SIZE, MIN_DEP, MAX_DEP, MIN_DEP_REG):
    bam_index = bam_in + ".bai"
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index)
    ## pipe
    for ctg in ctg_ls:
        bed_out1 = os.path.join(out_dir, ctg+"_clip.bed")   # path/chr1_clip.bed
        find_from_clip(bam_in, ctg, bed_out1, MIN_CLIP_NUM, MIN_CLIP_LEN)
        bed_out2 = os.path.join(out_dir, ctg+"_cov.bed")    # path/chr1_cov.bed
        # find_by_cov(bam_in, ctg, bed_out2, MIN_DEP, MIN_DEP_REG)

        ctg_dp_dir = dp_file_dir + "/" + ctg
        dp_file = os.path.join(dp_file_dir, ctg, ctg+".regions.bed.gz")
        if not os.path.isdir(ctg_dp_dir): os.makedirs(ctg_dp_dir)
        min_MQ = 20
        # min_len = 2000
        cmd = ["mosdepth", "-Q", str(min_MQ), "-b", str(DP_WIN_SIZE), "-c", ctg, ctg_dp_dir+"/"+ctg, bam_in]    # mosdepth -b 100 test/NC_19 ../step1_mapping/aln.sorted.BAM      # 指定winsize
        print("Running: {}".format(" ".join(cmd)))
        subprocess.check_call(" ".join(cmd), shell=True)
        dp_ls = get_dp_ls(dp_file)
        chr_len = bam_reader.get_reference_length(ctg)
        CHR_INFO = chr_info(ctg, chr_len)
        find_by_dp(dp_ls, DP_WIN_SIZE, CHR_INFO, MIN_DEP, MAX_DEP, bed_out2) # dp_ls, dp_win_size, CHR_INFO, MIN_DP, MAX_DP, bed_out
        

def run_find_candidate_parallel(bam_in, ctg_ls, out_dir, dp_file_dir, num_threads, MIN_CLIP_NUM, MIN_CLIP_LEN, DP_WIN_SIZE, MIN_DEP, MAX_DEP, MIN_DEP_REG):
    all_reference_ids = ctg_ls
    # all_reference_ids = [r for r in pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index).references]
    random.shuffle(all_reference_ids)
    chunk_size = len(all_reference_ids) // num_threads + 1
    threads = []
    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    for i in range(num_threads):
        contigs_list = all_reference_ids[i*chunk_size:(i+1)*chunk_size]
        if not contigs_list:
            continue
        threads.append(mp.Process(target=run_find_candidate, args=(bam_in, contigs_list, out_dir, dp_file_dir, MIN_CLIP_NUM, MIN_CLIP_LEN, DP_WIN_SIZE, MIN_DEP, MAX_DEP, MIN_DEP_REG)))

    signal.signal(signal.SIGINT, orig_sigint)
    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
            if t.exitcode != 0:
                raise Exception("One of the processes exited with code: {0}".format(t.exitcode))
    except KeyboardInterrupt:
        for t in threads:
            t.terminate()
        raise
    ## bed merge
    bed_ls = []
    for ctg in ctg_ls:
        bed_ls.append(os.path.join(out_dir, ctg+"_clip.bed"))
        bed_ls.append(os.path.join(out_dir, ctg+"_cov.bed"))
    bed_out_merge = os.path.join(out_dir, "candidate.bed")
    bed_merge(bed_ls, bed_out_merge)

def run_find_candidate2(bam_in, ctg, out_dir, config, dp_info_dic):
    # bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in + ".bai")
    ## pipe
    clip_params = config["clip_params"]
    min_clip_num = clip_params["min_clip_num"]
    min_clip_len = clip_params["min_clip_len"]
    bed_out1 = os.path.join(out_dir, "clip", ctg+"_clip.bed")   # path/chr1_clip.bed
    try:
        find_from_clip(bam_in, ctg, bed_out1, min_clip_num, min_clip_len)
    except Exception as e:
        print(f"Error: {e}")
    # 
    bed_out2 = os.path.join(out_dir, "cov", ctg+"_cov.bed")    # path/chr1_cov.bed
    dp_params = config["dp_params"]
    dp_info = dp_info_dic[ctg]
    try:
        find_by_dp2(dp_params, dp_info, bed_out2) # dp_ls, dp_win_size, CHR_INFO, MIN_DP, MAX_DP, bed_out
    except Exception as e:
        print(f"Error: {e}")
    print("{} find by dp clip Done!!!".format(ctg))

def write_dic(dic, fout):
    data=pd.DataFrame(dic)
    data.to_csv(fout, sep="\t", index=False)
def write_dpinfo_dic(dp_info_dic:dict, fout):
    data=pd.DataFrame(dp_info_dic)
    data.to_csv(fout, sep="\t", index=False)
    # with open(fout, "w") as f:
    #     f.write("chr_id\tchr_len\tchr_dp\twin_size\tblock_size\tblock_num\tblock_dp_ls\n")
    #     for dp_info in dp_info_dic.values():
    #         f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(dp_info.chr_id, dp_info.chr_len, dp_info.chr_dp, dp_info.win_size, dp_info.block_size, dp_info.block_num, dp_info.block_dp_ls, ))

def call_back(res):
    print(res)
def error_call_back(error_code):
    print(error_code)

def run_find_candidate2_parallel(ref, bam_in, ctg_ls, out_dir, threads, config):
    print("get_dp_info!!!")
    dp_win_size = config["dp_params"]["dp_win_size"]
    block_size = config["dp_params"]["block_size"]
    ## 
    dp_info_dic = get_dp_info_parallel(bam_in, threads, out_dir, dp_win_size, block_size)
    print("get_dp_info done!!!")

    ## 
    dirs = [os.path.join(out_dir, "clip"), os.path.join(out_dir, "cov")]
    for dir in dirs:
        if not os.path.exists(dir): os.makedirs(dir)
    pool = Pool(processes=threads)
    for ctg in ctg_ls:
        pool.apply_async(run_find_candidate2, args=(bam_in, ctg, out_dir, config, dp_info_dic), callback=call_back, error_callback=error_call_back) 
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    # for ctg in ctg_ls:
    #     run_find_candidate2(bam_in, ctg, out_dir, config, dp_info_dic)

    ## find from pileup
    if config["apply_pileup"]:
        print("apply pileup")
        pileup_dir = os.path.join(out_dir, "pileup")
        if not os.path.exists(pileup_dir):os.makedirs(pileup_dir)
        parse_pileup_parallel(ctg_ls, ref, bam_in, threads, config["pileup_params"], pileup_dir)
    
    print("find_candidate done!!!")


    ## bed merge
    bed_ls = []
    if config["apply_pileup"]:
        bed_ls.append(os.path.join(pileup_dir, "candidate_pileup.bed"))
    for ctg in ctg_ls:
        bed_ls.append(os.path.join(out_dir, "clip", ctg+"_clip.bed"))
        bed_ls.append(os.path.join(out_dir, "cov", ctg+"_cov.bed"))
    bed_out_merge = os.path.join(out_dir, "candidate.bed")
    # bed_merge(bed_ls, bed_out_merge)
    merge_cmd = ["cat", " ".join(bed_ls), "| sort -k 1,1 -k 2n,2", ">", bed_out_merge]
    subprocess.check_call(" ".join(merge_cmd), shell=True)


def main():
    ## Debug
    # bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/step1_mapping/aln.sorted.bam"
    # ctg_ls = ["NC_000019.10"]
    # out_dir = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/find_candidate_regions/test_find"
    # dp_file_dir = "/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/depths"
    # num_threads = 1
    # MIN_CLIP_NUM = 5
    # MIN_CLIP_LEN = 1000
    # DP_WIN_SIZE = 100
    # MIN_DEP = 10
    # MAXDEP = 48
    # MIN_DEP_REG = 100
    # run_find_candidate_parallel(bam_in, ctg_ls, out_dir, dp_file_dir, num_threads, MIN_CLIP_NUM, MIN_CLIP_LEN, DP_WIN_SIZE, MIN_DEP, MAXDEP, MIN_DEP_REG)
    return
    # ##
    # bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/resolve_region/solve_del/test/NC_060930.1_hap1.sort.bam"
    # bam_index = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/resolve_region/solve_del/test/NC_060930.1_hap1.sort.bam.bai"
    # # bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/resolve_region/solve_del/test/filter/filtered.bam"
    # # bam_index = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/resolve_region/solve_del/test/filter/filtered.bam.bai"
    # ctg_ls = ["NC_060930.1"]
    # bed_out1 = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/find_candidate_regions/test_find/find_from_clip.bed"

    # find_from_clip(bam_in, bam_index, ctg_ls, bed_out1)

    # ## 
    # fai = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/T2T_CHM13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.fai"
    # # depth = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/resolve_region/solve_del/test/find_from_filter/depths/NC_060930.1_hap1.depth"
    # # depth = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/resolve_region/solve_del/test/find_from_filter/depths2/NC_060930.1_hap1.depth"
    # depth = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/resolve_region/solve_del/test/test_depths/depths/NC_060930.1_hap1.depth"
    # ctg_ls = ["NC_060930.1"]
    # bed_out2 = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/find_candidate_regions/test_find/find_from_cov.bed"

    # find_by_cov(bam_in, bam_index, ctg_ls, bed_out2)


    # ## region_candidate
    # bed_out_merge = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/find_err/find_candidate_regions/find_candidate.bed"
    # bed_merge([bed_out1, bed_out2], bed_out_merge)


    ### ,MIN_CLIP_NUM, MIN_DEP, MIN_DEP_REG
    parser = argparse.ArgumentParser(description="get_candidate_regions")
    parser.add_argument("-t", "--threads", dest="threads", type=int, default=1)
    parser.add_argument("--ctg-ls", dest="ctg_ls", required=True, help="ctg list, format like:chr1,chr2")
    # parser.add_argument
    parser.add_argument("--bam", dest="bam_in", required=True, help="bam file")
    parser.add_argument("--out-dir", dest="out_dir", required=True, help="out dir of candidate regions")
    parser.add_argument("--min_clip_num", dest="min_clip_num", default=5, help="min_clip_num of window to be selected") # important
    parser.add_argument("--min_clip_len", dest="min_clip_len", default=1000, help="min_clip_len of clip to be selected as clip")
    parser.add_argument("--min-dep", dest="min_dep", default=3, help="dep threshold to be select as low dep region")
    parser.add_argument("--min-dep-reg", dest="min_dep_reg", default=100, help="minimum length of low dep region")
    parser.add_argument("-a", dest="all_chrs", action='store_true')
    args = parser.parse_args()
    if not os.path.exists(args.bam_in + ".bai"):
        print("BAM index is not available!!!")
        pysam.index("-@", str(args.threads), args.bam_in)
    if args.all_chrs:
        ctg_ls = list(pysam.AlignmentFile(args.bam_in, "rb").references)
    else:
        ctg_ls = args.ctg_ls.strip().split(",")        # 解析contig ls
    run_find_candidate_parallel(args.bam_in, ctg_ls, args.out_dir, int(args.threads), int(args.min_clip_num), int(args.min_clip_len), int(args.min_dep), int(args.min_dep_reg))
    
    # bed_ls = []
    # for ctg in ctg_ls:
    #     bed_ls.append(os.path.join(args.out_dir, ctg+"_clip.bed"))
    #     bed_ls.append(os.path.join(args.out_dir, ctg+"_cov.bed"))
    # bed_out_merge = os.path.join(args.out_dir, "candidate.bed")
    # bed_merge(bed_ls, bed_out_merge)

    # bam_in = sys.argv[1]
    # bam_index = sys.argv[2]
    # out_dir = sys.argv[3]
    # ctg_ls = sys.argv[4].strip().split(",")        # 解析contig ls
    # num_threads = int(sys.argv[5])
    ## single pipe
    # run_find_candidate(bam_in, bam_index, ctg_ls, out_dir)
    ## parallel
    # run_find_candidate_parallel(bam_in, bam_index, ctg_ls, out_dir, num_threads)

    ### merge bed
    # bed_ls = []
    # for ctg in ctg_ls:
    #     bed_ls.append(os.path.join(out_dir, ctg+"_clip.bed"))
    #     bed_ls.append(os.path.join(out_dir, ctg+"_cov.bed"))
    # bed_out_merge = os.path.join(out_dir, "candidate.bed")
    # bed_merge(bed_ls, bed_out_merge)
    print("------------------------merge candidate bed Done--------------------------")

if __name__ == "__main__":
    print("Begin")
    import yaml
    bam_in = "/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe/step1_mapping/aln.sorted.bam"
    ctg_ls = ["NC_000018.10"]
    out_dir = "/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe/test2"
    threads = 2
    config_file = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Configs/Config.yaml"
    with open(config_file, "r")as f:
        config = yaml.safe_load(f.read())
    print(config)
    run_find_candidate2_parallel(bam_in, ctg_ls, out_dir, threads, config)
    pass
