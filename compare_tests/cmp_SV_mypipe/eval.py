


# import pandas as pd
# import math
from collections import defaultdict
'''
asm_patch
asm_patch,N_fill,asm_patch
asm_patch,N_fill,N_fill
N_fill
N_fill,N_fill,asm_patch
reads_patch
'''
'''
reads_patch
asm_patch
N_fill
'''

# 打开区间数据文件并读取数据
# @DeprecationWarning
def parse_bed_info(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # 创建一个字典来存储区间信息
    interval_data1 = {}

    for line in lines:
        if line.startswith("#"): continue
        parts = line.strip().split('\t')  # 假设数据以制表符分隔
        # if len(parts) == 3:
        chr_id, start, end = parts[:3]
        start, end = int(start), int(end)

        # 存储区间信息
        if chr_id not in interval_data1:
            interval_data1[chr_id] = []
        interval_data1[chr_id].append((start, end))
    
    for intervals in interval_data1.values():
        intervals.sort(key=lambda reg:reg[0])
    # print(interval_data1["NC_000001.11"])
    return interval_data1


def parse_consensus_bed_info(file_path, keys):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # 创建一个字典来存储区间信息
    # interval_data0 = dict.fromkeys(keys, [])
    # interval_data1 = dict.fromkeys(keys, [])
    # interval_data2 = dict.fromkeys(keys, [])
    # interval_data3 = dict.fromkeys(keys, [])
    interval_data0 = defaultdict(list) # Sum
    interval_data1 = defaultdict(list) # reads_patch
    interval_data2 = defaultdict(list) # asm_patch
    interval_data3 = defaultdict(list) # N_fill

    for line in lines:
        if line.startswith("#"): continue
        parts = line.strip().split('\t')  # 假设数据以制表符分隔
        # if len(parts) == 3:
        chr_id, start, end, ops, info, patch_id = parts[:6]
        start, end = int(start), int(end)
        op_ls = ops.split(",")
        # 存储区间信息
        interval_data0[chr_id].append((start, end)) # all

        if len(op_ls) > 1:  # N_fill + asm_patch
            # mid = (end - start) // 2 + start
            # interval_data2[chr_id].append((start, mid))
            # interval_data3[chr_id].append((mid, end))
            interval_data2[chr_id].append((start, end))
        else:   # asn_patch and 
            if op_ls[0] == "reads_patch":
                interval_data1[chr_id].append((start, end))
            elif op_ls[0] == "asm_patch":
                interval_data2[chr_id].append((start, end))
            elif op_ls[0] == "N_fill":
                interval_data3[chr_id].append((start, end))
            else: raise ValueError
    # print(interval_data1["NC_000001.11"])
    return interval_data0, interval_data1, interval_data2, interval_data3


def cmp_chr_reg():
    
    pass

def stats_chr_reg(intervals):   # 输入区间数组，进行数据统计：区间数目、长度
    # 计算单组的统计信息
    num_intervals = len(intervals)
    total_length = sum(abs(end - start) for start, end in intervals)
    average_length = total_length / num_intervals
    longest_interval = max(intervals, key=lambda x: x[1] - x[0])
    shortest_interval = min(intervals, key=lambda x: x[1] - x[0])
    
    return {
        # 'num_intervals': num_intervals,
        'total_length': total_length // 1000,
        # 'average_length': average_length,
        # 'longest_interval': longest_interval,
        # 'shortest_interval': shortest_interval
    }

def stats_all(interval_dic:defaultdict):    # {chr:[[], []], }
    ls = []
    print("Stats All")
    # print(interval_dic)
    for chr, intervals in interval_dic.items():
        stats_dic = stats_chr_reg(intervals)
        ls.append(chr)
        Print_KV(stats_dic)
    print(ls)

def Print_KV(stats_dic):
    # for key, values in stats_dic.items():
    #     print("{}:\t{}".format(key, values))
    # for key in stats_dic.keys():
    #     print(key)
    for val in stats_dic.values():
        print(val)

print("**********SV tools***********")
file_path = "/public/home/hpc214712170/Test/tests/chm13_ont/my_pipe2/SV/cute_SV/For_cmp/merge.out.bed"  # 替换为实际文件路径
interval_data = parse_bed_info(file_path)
stats_dic = stats_chr_reg(interval_data["NC_000001.11"])
Print_KV(stats_dic)

# 
print("**********Mypipe***********")
keys = ['NC_000001.11', 'NC_000002.12', 'NC_000003.12', 'NC_000004.12', 'NC_000005.10', 'NC_000006.12', 'NC_000007.14', 'NC_000008.11', 'NC_000009.12', 'NC_000010.11', 'NC_000011.10', 'NC_000012.12', 'NC_000013.11', 'NC_000014.9', 'NC_000015.10', 'NC_000016.10', 'NC_000017.11', 'NC_000018.10', 'NC_000019.10', 'NC_000020.11', 'NC_000021.9', 'NC_000022.11', 'NC_000023.11']
consensus_file_path = "/public/home/hpc214712170/Test/tests/chm13_ont/my_pipe2/step3_SV_consensus/consensus.bed"
interval_data0, interval_data1, interval_data2, interval_data3 = parse_consensus_bed_info(consensus_file_path, keys)

stats_dic0 = stats_chr_reg(interval_data0["NC_000001.11"]) 
stats_dic1 = stats_chr_reg(interval_data1["NC_000001.11"])
stats_dic2 = stats_chr_reg(interval_data2["NC_000001.11"])
stats_dic3 = stats_chr_reg(interval_data3["NC_000001.11"])
# print(stats_dic0)

Print_KV(stats_dic0)
Print_KV(stats_dic1)
Print_KV(stats_dic2)
Print_KV(stats_dic3)


# stats_all(interval_data0)
stats_all(interval_data1)
stats_all(interval_data2)
stats_all(interval_data3)

