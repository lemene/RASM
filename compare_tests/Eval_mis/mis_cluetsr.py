
import sys
from collections import defaultdict

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
    bed_in = sys.argv[1]
    bed_out = sys.argv[2]
    cluster_len = int(sys.argv[3])  # 1000
    dic = defaultdict(list)
    new_dic = defaultdict(list)
    with open(bed_in, "r") as f:
        for line in f:
            if line.startswith("#"): continue
            if not line: continue
            ls = line.strip().split()
            if ";" in ls[1]:
                # continue
                ctg, start_, end_ = ls[:3]
                start_ls = start_.split(";")
                end_ls = end_.split(";")
                start_ls = [int(i) for i in start_ls]
                end_ls = [int(i) for i in end_ls]
                start = min(start_ls)
                end = max(end_ls)
            else:
                ctg, start, end = ls[:3]
                start, end = int(start), int(end)
            dic[ctg].append([start, end])
    mis_num = 0
    raw_mis_num = 0
    for ctg, mis_ls in dic.items():
        raw_mis_num += len(mis_ls)
        cluster_ls = cluster(mis_ls, cluster_len)
        new_dic[ctg] = cluster_ls
        mis_num += len(cluster_ls)
    with open(bed_out, "w") as f:
        f.write("#ctg\tstart\tend\n")
        f.write("#cluster_dis={}\n".format(cluster_len))
        f.write("#raw_mis_num={}\n".format(raw_mis_num))
        f.write("#mis_num={}\n".format(mis_num))
        for ctg, ls in new_dic.items():
            for reg in ls:
                f.write("{}\t{}\t{}\n".format(ctg, str(reg[0]), str(reg[1])))

# python /public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/mis_cluetsr.py mis.bed mis_cl.bed 5000