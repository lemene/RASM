
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

class Depth_info():

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