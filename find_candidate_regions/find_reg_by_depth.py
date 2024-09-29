import os
import numpy as np
import gzip
from multiprocessing import Pool
import subprocess
import pysam
import time
from collections import namedtuple
# from io import TextIOWrapper
# import matplotlib.pyplot as plt
chr_info = namedtuple('chr_info', ["chr_id", "chr_len"])
Region = namedtuple('Region', ["chr_id", "start", "end"])

######################## Depth parse ##############################
class Depth_info():
    '''
    cal some depth stastic of depth
    计算覆盖度的统计数据 per contig
    '''
    def __init__(self, chr_id, chr_len, dp_ls, dp_win_size, block_size, whole_dp) -> None:   # win_size表示输入dp_ls采用的dp
        # change_ratio = 0.5  # 波动范围在均值的这个之内的算作是正常的覆盖度
        min_cover_ratio = 0.5   # 正常覆盖度值的区间比例>这个值，适用局部计算的depth，否则，适用全局或周边的组为depth
        
        self.dp_ls = np.array(dp_ls)
        self.win_size = dp_win_size
        self.block_size = block_size
        self.chr_len = chr_len
        self.chr_id = chr_id
        self.whole_dp = whole_dp
        ## 
        block_batch = self.block_size//self.win_size # block包含多少个dp_ls的window
        self.block_dp_ls = []     # 
        self.chr_dp = whole_dp  # 初始化 
        self.block_num = self.chr_len // self.block_size + 1 if chr_len % block_size > 0 else chr_len // block_size  # 有多少个block

        ### 排除异常元素排除了0，由于基因组有N填充的区域，保留
        # criteria = (self.dp_ls > 0.01) & (self.dp_ls < 60)  
        criteria = self.dp_ls > 0.01    # 保留元素的条件, 避免异常值的影响（尤其是0值的影响）
        new_data = np.extract(criteria, self.dp_ls)
        if self.dp_ls.size > 0:
            self.cov_ratio = new_data.size / self.dp_ls.size
        else:
            self.cov_ratio = -1
        self.chr_avg_dp = np.average(new_data) 
        if len(new_data) / len(self.dp_ls) < 0.8 or np.median(new_data) < whole_dp * 0.9 or np.median(new_data) > whole_dp * 1.1 or self.block_num < 2:   # 当前contig计算的depth不可信，使用所有contig的depth均值来替代
            print("Error contig depth: {}".format(self.chr_id))
            self.chr_dp = whole_dp
            self.block_dp_ls = [whole_dp] * self.block_num
        else:
            low_dp, high_dp = np.quantile(new_data, [0.10, 0.90])   # dp to select normalr dp to calculate window avg dp
            self.chr_dp = np.median(new_data)   # 全局的dp
            print("{} dp Range:{}-{}".format(self.chr_id, low_dp, high_dp))
            ## compute dp of a block
            for i in range(self.block_num):
                ls = [] # ls 统计一个block内符合条件的
                if self.block_size * (i + 1) < chr_len:   # block: 
                    for j in range(block_batch * i, block_batch * (i + 1)):
                        if self.dp_ls[j] < high_dp and self.dp_ls[j] > low_dp:
                            ls.append(self.dp_ls[j])
                    if len(ls) / block_batch > min_cover_ratio:
                        self.block_dp_ls.append(np.median(ls))
                    else:
                        print("{}-{}: use whole contig dp val instead".format(self.block_size * i, self.block_size * (i + 1)))
                        self.block_dp_ls.append(self.chr_dp)   # 用全局的代表
                else:
                    self.block_dp_ls.append(self.chr_dp)  # 适用全局的代表把
    @staticmethod
    def write_dp_info(dpinfo_dic:dict, fout):
        head_ls = ["chr_id", "chr_len", "whole_dp", "chr_avg_dp", "chr_dp", "cov_ratio", \
                   "win_size", "block_size", "block_num", "block_dp_ls"]
        with open(fout, "w") as f:
            f.write("#" + "\t".join(head_ls) + "\n")
            for dp_info in dpinfo_dic.values():
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
                dp_info.chr_id, dp_info.chr_len, dp_info.whole_dp, dp_info.chr_avg_dp, \
                dp_info.chr_dp, dp_info.cov_ratio, dp_info.win_size, dp_info.block_size, \
                dp_info.block_num, ",".join(str(dp) for dp in dp_info.block_dp_ls)))
    @staticmethod
    def read_dp_info(file):
        head_ls = ["chr_id", "chr_len", "whole_dp", "chr_avg_dp", "chr_dp", "cov_ratio", \
                   "win_size", "block_size", "block_num", "block_dp_ls"]
        data_ls = []
        dp_info_dic = {}
        with open(file, "r") as f:
            for line in f:
                ctg_info_dic = {}
                if not line: continue 
                if line.startswith("#"):
                    continue
                data_ls = line.strip().split("\t")
                for i in range(len(head_ls)):
                    ctg_info_dic[head_ls[i]] = data_ls[i]
                ctg_info_dic["chr_len"] = int(ctg_info_dic["chr_len"])
                ctg_info_dic["whole_dp"] = float(ctg_info_dic["whole_dp"])
                ctg_info_dic["chr_avg_dp"] = float(ctg_info_dic["chr_avg_dp"])
                ctg_info_dic["chr_dp"] = float(ctg_info_dic["chr_dp"])
                ctg_info_dic["cov_ratio"] = float(ctg_info_dic["cov_ratio"])
                ctg_info_dic["win_size"] = int(ctg_info_dic["win_size"])
                ctg_info_dic["block_size"] = int(ctg_info_dic["block_size"])
                ctg_info_dic["block_num"] = int(ctg_info_dic["block_num"])
                ctg_info_dic["block_dp_ls"] = [float(dp) for dp in ctg_info_dic["block_dp_ls"].split(",")]
                ## 
                dp_info_dic[ctg_info_dic["chr_id"]] = ctg_info_dic
        return dp_info_dic

def get_Depth_info(chr_id, chr_len, dp_ls, dp_win_size, block_size, whole_dp):
    return Depth_info(chr_id, chr_len, dp_ls, dp_win_size, block_size, whole_dp)


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

def run_mosdepth(ctg, bam_in, out_dir, DP_WIN_SIZE):
    print("Run mosdepth for {}".format(ctg))
    min_MQ = 20
    # DP_WIN_SIZE = 100
    prefix = out_dir+"/"+ctg
    cmd = ["mosdepth", "-Q", str(min_MQ), "-b", str(DP_WIN_SIZE), "-c", ctg, prefix, bam_in]    # mosdepth -b 100 test/NC_19 ../step1_mapping/aln.sorted.BAM      # 指定winsize
    print("Running: %s", " ".join(cmd))
    subprocess.check_call(" ".join(cmd), shell=True)
    print("Run {} done".format(ctg))

def get_dp(ctg, bam_in, ctg_out_dir, DP_WIN_SIZE, min_MQ):
    '''per ctg'''
    if not os.path.isdir(ctg_out_dir):os.makedirs(ctg_out_dir)
    dp_out_file = ctg_out_dir+"/"+ctg + ".regions.bed.gz"
    ## 
    run_mosdepth2(ctg, bam_in, ctg_out_dir, DP_WIN_SIZE, 4, min_MQ)
    dp_ls = read_mosdepth_dp_file(dp_out_file)
    return {ctg: dp_ls}     # 返回一个字典存储

def get_dp_info_parallel(bam_in, threads, out_dir, DP_WIN_SIZE, Block_size, min_MQ):
    print("----------------get_dp_info_parallel----------------")
    dp_file_dir = os.path.join(out_dir, "depths")
    dp_info_dir = os.path.join(out_dir, "dp_info")
    if not os.path.isdir(dp_file_dir):os.makedirs(dp_file_dir)
    if not os.path.isdir(dp_info_dir):os.makedirs(dp_info_dir)

    bam_reader = pysam.AlignmentFile(bam_in, "rb") 
    ctg_ls = bam_reader.references
    dp_dic = {}
    pool = Pool(processes=threads)
    results = [pool.apply_async(get_dp, args=(ctg, bam_in, dp_file_dir + "/" + ctg, DP_WIN_SIZE, min_MQ)) for ctg in ctg_ls]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕

    for res in results:
        dp_dic.update(res.get())    # 更新所有的

    whole_dp_ls = []
    for ctg, ctg_dp_ls in dp_dic.items():
        whole_dp_ls.extend(ctg_dp_ls)
    whole_dp = np.median(whole_dp_ls)   # 全局的dp，使用中位数表示
    print("Whole dp: {} !!!".format(whole_dp))

    ## 
    dpinfo_dic = {}
    
    results = []
    pool = Pool(processes=threads)
    for ctg, ctg_dp_ls in dp_dic.items():
        ctg_len = bam_reader.get_reference_length(ctg)
        results.append(pool.apply_async(get_Depth_info, args=(ctg, ctg_len, ctg_dp_ls, DP_WIN_SIZE, Block_size, whole_dp)))
        # results = [pool.apply_async(get_Depth_info, args=(ctg, ctg_len, ctg_dp_ls, DP_WIN_SIZE, Block_size, whole_dp))]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    for res in results:
        ctg_dpinfo = res.get()
        dpinfo_dic[ctg_dpinfo.chr_id] = ctg_dpinfo
    
    dp_info_file = dp_info_dir + "/" + "dpinfo.bed"
    Depth_info.write_dp_info(dpinfo_dic, dp_info_file)
    return dpinfo_dic


####################################################################
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

def find_by_dp2(dp_params, dp_info:Depth_info, bed_out):  # 传如dp列表，以及dp计算的window大小
    '''
    find depth reg by cov, 
    
    '''
    dp_lower_bound = dp_params["dp_lower_bound"]
    dp_upper_bound = dp_params["dp_upper_bound"]
    chr_id = dp_info.chr_id
    chr_len = dp_info.chr_len
    dp_ls = dp_info.dp_ls
    dp_win_size = dp_info.win_size
    block_size = dp_info.block_size
    block_num = dp_info.block_num
    block_batch = block_size // dp_win_size # 每个block有多少个window
    ## 对每个block使用block_dp作为标准值，根据每个window的dp值，计算偏移值
    ls = []
    for i in range(block_num):  # 块范围：[i*block_size, (i+1)*block_size]
        block_dp_ls = dp_info.block_dp_ls
        block_upper_dp = block_dp_ls[i] * dp_upper_bound
        block_lower_dp = block_dp_ls[i] * dp_lower_bound
        print("Block {}:{}-{}, dp limit:{}-{}".format(chr_id, i * block_size, (i+1)*block_size, block_lower_dp, block_upper_dp))
        for j in range(block_batch):    # window范围：[]
            if i*block_size + (j + 1)*dp_win_size < chr_len:
                win_l = i*block_size + j*dp_win_size
                win_r = i*block_size + (j + 1)*dp_win_size
                win_dp = dp_ls[i*block_batch + j]   # 窗口dp
                if win_dp < block_lower_dp or win_dp > block_upper_dp:
                    ls.append(Region(chr_id, win_l, win_r))
            elif i*block_size + j*dp_win_size < chr_len:    # 最后一个window
                win_l = i*block_size + j*dp_win_size
                win_r = chr_len
                win_dp = dp_ls[i*block_batch + j]   # 窗口dp
                if win_dp < block_lower_dp or win_dp > block_upper_dp:
                    ls.append(Region(chr_id, win_l, win_r))
                break
            else:   # i*block_size + j*dp_win_size >= chr_len
                break
    cluster_dis = dp_params["cluster_dis"]
    final_ls = cluster_by_dis(ls, cluster_dis)

    with open(bed_out, "w") as fout:
        for reg in final_ls:
            fout.write("{}\t{}\t{}\t{}bp-cov_reg\n".format(reg.chr_id, reg.start, reg.end, reg.end-reg.start))
            # print("{}\t{}\t{} cov_reg".format(reg.chr_id, reg.start, reg.end))
    return final_ls