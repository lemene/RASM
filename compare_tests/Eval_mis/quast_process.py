
'''
1、排除quast中十分密集的Misassembly干扰
2、排除小的组装碎片的干扰
'''

import sys
# from Eval__mis_SV.convert2 import Mis_rec
import convert2
def read_fai(fai_filename):
    chromosome_lengths = {}
    try:
        with open(fai_filename, 'r') as fai_file:
            for line in fai_file:
                parts = line.strip().split('\t')  # Split each line into parts using tab as a delimiter
                if len(parts) >= 2:
                    chrom_name = parts[0]  # The first part is the chromosome name
                    chrom_length = int(parts[1])  # The second part is the chromosome length
                    chromosome_lengths[chrom_name] = chrom_length
    except IOError:
        print(f"Error: File {fai_filename} does not exist or cannot be read.")
    return chromosome_lengths


def cluster(intervals, dis):
    if not intervals:
        return []
    # 初始化结果列表，以及当前合并区间的起始和结束位置
    merged = []
    current_start, current_end = intervals[0]

    for i in range(1, len(intervals)):
        interval_start, interval_end = intervals[i]
        
        # 当前区间的起始位置与上一个区间的结束位置的距离小于等于指定距离时，进行合并
        if interval_start - current_end <= dis:
            # 合并区间的新结束位置为当前区间和上一个合并区间中的较大值
            current_end = max(current_end, interval_end)
        else:
            # 将上一个合并后的区间添加到结果列表中
            merged.append((current_start, current_end))
            # 重置合并区间的起始和结束位置为当前区间的值
            current_start, current_end = interval_start, interval_end

    # 添加最后一个合并后的区间到结果列表中
    merged.append((current_start, current_end))

    return merged

if __name__ == "__main__":

    fai = sys.argv[1]
    bed_in = sys.argv[2]
    bed_out = sys.argv[3]
    min_ctg = int(sys.argv[4])   # 1000000
    cluster_len = int(sys.argv[5])  # 10000

    Mis_rec_dic = convert2.Mis_rec.read_rec(bed_in)
    ctg_len_dic = read_fai(fai)

    new_dic = {}
    skip_ls = []
    for ctg, ctg_len in ctg_len_dic.items():
        if ctg_len_dic[ctg] < min_ctg:
            print("Skip:{}".format(ctg))
            skip_ls.append(ctg)
    for ctg, ctg_rec_ls in Mis_rec_dic.items():
        if ctg in skip_ls: continue
        # ctg_rec_ls.sort(key=lambda rec: rec.start)
        ctg_rec_ls.sort()
        cluster_ls = cluster(ctg_rec_ls, cluster_len)
        print("Cluster form {}\n->{}".format(ctg_rec_ls, cluster_ls))
        new_dic[ctg] = cluster_ls

    f = open(bed_out, "w")
    for ctg, ls in new_dic.items():
        for reg in ls:
            f.write("{}\t{}\t{}\n".format(ctg, str(reg[0]), str(reg[1])))
    f.close()

    with open(bed_out + ".ex", "w") as f:
        for ctg in skip_ls:
            f.write("{}\n".format(ctg))