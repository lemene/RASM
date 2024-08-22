# import os
# import sys
import numpy as np
import gzip
from multiprocessing import Pool
import subprocess
import os
import pysam
import time

def run_mosdepth2(out_dir, bam_in, threads):
    min_MQ = 20
    DP_WIN_SIZE = 100
    cmd = ["mosdepth", "-Q", str(min_MQ), "-b", str(DP_WIN_SIZE), "-t", str(threads), out_dir + "/all/" + "all", bam_in]    # mosdepth -b 100 test/NC_19 ../step1_mapping/aln.sorted.BAM      # 指定winsize
    print("Running: %s", " ".join(cmd))
    subprocess.check_call(" ".join(cmd), shell=True)

def read_mosdepth_dp_file(dp_file):
    '''best for ctg dp file'''
    dp_ls = []
    with gzip.open(dp_file, "rt") as f:
        for line in f:
            fields = line.strip().split("\t")
            # dp = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
            dp = float(fields[3])
            dp_ls.append(dp)
    return dp_ls

class Depth_info():
    '''
    cal some depth stastic of depth
    计算覆盖度的统计数据
    '''
    def __init__(self, chr_id, chr_len, dp_ls, dp_win_size, block_size, whole_dp) -> None:   # win_size表示输入dp_ls采用的dp
        # change_ratio = 0.5  # 波动范围在均值的这个之内的算作是正常的覆盖度
        min_cover_ratio = 0.5   # 正常覆盖度值的区间比例>这个值，适用局部计算的depth，否则，适用全局或周边的组为depth
        
        self.dp_ls = np.array(dp_ls)
        self.win_size = dp_win_size
        self.block_size = block_size
        self.chr_len = chr_len
        self.chr_id = chr_id
        ## 
        block_batch = self.block_size//self.win_size # block包含多少个dp_ls的window
        self.block_dp = []     # 
        self.chr_dp = whole_dp  # 初始化
        self.block_num = self.chr_len // self.block_size + 1 if chr_len % block_size > 0 else chr_len // block_size  # 有多少个block

        ### 排除异常元素排除了0，由于基因组有N填充的区域，保留
        # criteria = (self.dp_ls > 0.01) & (self.dp_ls < 60)  
        criteria = self.dp_ls > 0.01    # 保留元素的条件, 避免异常值的影响（尤其是0值的影响）
        new_data = np.extract(criteria, self.dp_ls)
        if len(new_data) / len(self.dp_ls) < 0.5:   # 档期那contig计算的depth不可信，使用所有contig的depth均值来替代
            print("Error contig depth: {}".format(self.chr_id))
            self.chr_dp = whole_dp
            self.block_dp = [whole_dp] * self.block_num
        else:
            low_dp, high_dp = np.quantile(new_data, [0.10, 0.90])   # dp to select normalr dp to calculate window avg dp
            self.chr_dp = np.median(new_data)   # 全局的dp
            print("Range:{}-{}".format(low_dp, high_dp))
            ## compute dp of a block
            for i in range(self.block_num):
                ls = [] # ls 统计一个block内符合条件的
                if self.block_size * (i + 1) < chr_len:   # block: 
                    for j in range(block_batch * i, block_batch * (i + 1)):
                        if self.dp_ls[j] < high_dp and self.dp_ls[j] > low_dp:
                            ls.append(self.dp_ls[j])
                    if len(ls) / block_batch > min_cover_ratio:
                        self.block_dp.append(np.median(ls))
                    else:
                        print("{}-{}: use whole contig dp val instead".format(self.block_size * i, self.block_size * (i + 1)))
                        self.block_dp.append(self.chr_dp)   # 用全局的代表
                else:
                    self.block_dp.append(self.chr_dp)  # 适用全局的代表把

# def cal_depth_info():
def run_mosdepth(ctg, bam_in, out_dir, DP_WIN_SIZE):
    print("Run mosdepth for {}".format(ctg))
    min_MQ = 20
    # DP_WIN_SIZE = 100
    ctg_dp_dir = out_dir + "/" + ctg
    if not os.path.isdir(ctg_dp_dir):os.makedirs(ctg_dp_dir)
    prefix = out_dir+"/"+ctg
    cmd = ["mosdepth", "-Q", str(min_MQ), "-b", str(DP_WIN_SIZE), "-c", ctg, prefix, bam_in]    # mosdepth -b 100 test/NC_19 ../step1_mapping/aln.sorted.BAM      # 指定winsize
    print("Running: %s", " ".join(cmd))
    subprocess.check_call(" ".join(cmd), shell=True)
    print("Run {} done".format(ctg))

def get_dp_info(ctg, bam_in, ctg_out_dir, DP_WIN_SIZE):
    if not os.path.isdir(ctg_out_dir):os.makedirs(ctg_out_dir)
    dp_out_file = ctg_out_dir+"/"+ctg + ".regions.bed.gz"
    ## 
    run_mosdepth(ctg, bam_in, ctg_out_dir, DP_WIN_SIZE)
    dp_ls = read_mosdepth_dp_file(dp_out_file)
    return {ctg: dp_ls}     # 返回一个字典存储

def get_dp_info_parallel(bam_in, threads, out_dir, DP_WIN_SIZE, Block_size):
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai") 
    ctg_ls = bam_reader.references
    dp_dic = {}
    pool = Pool(processes=threads)
    results = [pool.apply_async(get_dp_info, args=(ctg, bam_in, out_dir + "/" + ctg, DP_WIN_SIZE)) for ctg in ctg_ls]
    pool.close() # 关闭进程池，表示不能再往进程池中添加进程，需要在join之前调用
    pool.join() # 等待进程池中的所有进程执行完毕
    for res in results:
        dp_dic.update(res.get())    # 更新所有的
    print("Get dp ls done")
    whole_dp_ls = []
    for ctg, ctg_dp_ls in dp_dic.items():
        whole_dp_ls.extend(ctg_dp_ls)
    whole_dp = np.median(whole_dp_ls)   # 全局的dp
    print("Whole dp: {} !!!".format(whole_dp))
    ## 
    dp_info_dic = {}
    for ctg, ctg_dp_ls in dp_dic.items():
        ctg_len = bam_reader.get_reference_length(ctg)
        dp_info_dic[ctg] = Depth_info(ctg, ctg_len, ctg_dp_ls, DP_WIN_SIZE, Block_size, whole_dp)
    return dp_info_dic

if __name__ == "__main__":
    # # ls = ["A", "B"]
    # # print(",".join(ls))
    # import gzip
    # dp_file = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_19/my_pipe/depths/NC_000019.10/NC_000019.10.regions.bed.gz"
    # # dp_file = "/public/home/hpc214712170/Test/tests/all_data/chm13/hifi/depths/all/"
    # chr_len = 0
    # dp_ls = []
    # with gzip.open(dp_file, "rt") as f:
    #     for line in f:
    #         dp = float(line.strip().split()[3])
    #         dp_ls.append(dp)
    # dp_info = Depth(dp_ls, 100, 58617616, 5000000)
    # print("chr_dp:", dp_info.chr_dp)
    # print(dp_info.block_dp)
    # for i, dp in enumerate(dp_info.block_dp):
    #     if dp > dp_info.chr_dp * 1.3 or dp < dp_info.chr_dp * 0.7:
    #         print(dp, i*100000, (i+1)*100000)

    
    ## 
    t1 = time.time()
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe/step1_mapping/aln.sorted.bam"
    out_dir = "/public/home/hpc214712170/test/depths2"
    # threads = 24
    DP_WIN_SIZE = 100
    Block_size = 5000000
    # dic = get_dp_info_parallel(bam_in, threads, out_dir, DP_WIN_SIZE, Block_size)

    dp_dic = {}
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai") 
    ctg_ls = bam_reader.references
    # ctg_ls = ["NC_000022.11"]
    for ctg in ctg_ls:
        ctg_dp_ls = read_mosdepth_dp_file(out_dir + "/" + ctg + "/" + ctg + ".regions.bed.gz")
        dp_dic[ctg] = ctg_dp_ls
    t2 = time.time()
    print("Read file cost {}".format(t2 - t1))
    whole_dp_ls = []
    for ctg, ctg_dp_ls in dp_dic.items():
        whole_dp_ls.extend(ctg_dp_ls)
    whole_dp = np.median(whole_dp_ls)   # 全局的dp

    print("whole_dp:", whole_dp)
    whole_dp_ls = np.array(whole_dp_ls)
    criteria = whole_dp_ls > 0.01
    new_whole_dp_ls = np.extract(criteria, whole_dp_ls)
    whole_dp = np.median(new_whole_dp_ls)   # 全局的dp
    print("whole_dp:", whole_dp)
    ## 
    t3 = time.time()
    print("Cost:", t3 - t2)
    dp_info_dic = {}
    for ctg, ctg_dp_ls in dp_dic.items():
        ctg_len = bam_reader.get_reference_length(ctg)
        dp_info_dic[ctg] = Depth_info(ctg_dp_ls, DP_WIN_SIZE, ctg_len, Block_size, whole_dp)
    t4 = time.time()
    print("Cost {}".format(t4 - t3))
    pass
