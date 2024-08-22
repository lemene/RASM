import pysam
from collections import namedtuple
Region = namedtuple('Region', ["chr_id", "start", "end"])
def cal_sum(arry, l, r):
    return sum(arry[l:r])

def cal_avg_depth2(bam_in, reg):  # cal depth of a contig, return list of depth, length == contig length
    reg = reg.chr_id + ":" + str(reg.start) + "-" + str(reg.end)
    '''time samtools depth -a -Q $min_MQ -r $REG ${bam_file} > $out_prefix.depth    
    # 不加-J del处depth为0'''
    min_MQ = 40
    # reg_len = int(reg.split(":")[2]) - int(reg.split(":")[1])
    depth_ls = []
    depth_stream = pysam.depth("-a", "-Q", str(min_MQ), "-@", "5", "-r", reg, bam_in)    # depth_stream接收stdout
    for line in depth_stream.split("\n"):
        line_ls = line.split()
        if not line_ls:
            continue
        depth_ls.append(int(line_ls[2]))
    return cal_sum(depth_ls, 0, len(depth_ls)) / len(depth_ls)
def test_cal():
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/yeast/test_snakemake/snakemake_res/CRR198362_40X/step1_mapping/aln.sort.bam"
#     NC_001134.8     258000  260500  clip_reg
# NC_001134.8     265500  267000
# chrXI:75430-79067
    reg = Region(chr_id='chrXI', start=75430, end=79067)
    print(cal_avg_depth2(bam_in, reg))
    return

test_cal()

def is_INV(reg1, reg2, bam_in):

    return False

def is_pair(reg1, reg2, bam_in):    # 一前一后两个区间，判断是否为一对pair。根据bam来判断。main for dels
    if reg1.chr_id!= reg2.chr_id:
        return False
    max_connect_len = 20000      # 20000
    if reg1.start - reg2.end > max_connect_len:
        return False
    ## 根据区间内的信息判断: 1、读数支持：左右两边有相同的读数支持 2、depth支持：将删除不算作depth，计算区间depth
    MIN_SUPPOPRT_NUM = 10   # 10
    MIN_MAPPING_QUALITY = 20
    MIN_DP = 3   # 删除的深度
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in + ".bai")
    reads_set1 = set()
    reads_set2 = set()
    # print("-------------", reg1, reg2, "-------------")
    for read in bam_reader.fetch(reg1.chr_id, reg1.start, reg1.end):
        # print(read.query_name)
        if read.is_secondary and read.mapping_quality < MIN_MAPPING_QUALITY:
            continue
        reads_set1.add(read.query_name)
    for read in bam_reader.fetch(reg2.chr_id, reg2.start, reg2.end):
        if read.is_secondary and read.mapping_quality < MIN_MAPPING_QUALITY:
            continue
        reads_set2.add(read.query_name)
    commons = reads_set1.intersection(reads_set2)
    # print(commons)
    if cal_avg_depth2(bam_in, Region(reg1.chr_id, reg1.end, reg2.start)) < MIN_DP: # reg1.chr_id+":"+str(reg1.end)+"-"+str(reg2.start)
            print("pairs:{},{}-{} and {}-{}; support_reads:{}".format(reg1.chr_id, reg1.start, reg1.end, reg2.start, reg2.end, len(commons)))
            return True
    # if len(commons) >= MIN_SUPPOPRT_NUM:
        # if cal_avg_depth2(bam_in, Region(reg1.chr_id, reg1.end, reg2.start)) < MIN_DP: # reg1.chr_id+":"+str(reg1.end)+"-"+str(reg2.start)
        #     print("pairs:{},{}-{} and {}-{}; support_reads:{}".format(reg1.chr_id, reg1.start, reg1.end, reg2.start, reg2.end, len(commons)))
        #     return True
        # print("pairs:{},{}-{} and {}-{}; support_readpairs:{}".format(reg1.chr_id, reg1.start, reg1.end, reg2.start, reg2.end, len(commons)))
        # return True
    return False

# NC_001134.8     258000  259500  clip_reg
# NC_001134.8     265500  267000  clip_reg
# NC_001139.9     730000  731500  clip_reg
# NC_001139.9     735000  737000  clip_reg
# reg1 = Region("NC_001134.8", 258500, 259500)
# reg2 = Region("NC_001134.8", 265500, 267500)  # DEL   can be detected by depth
# reg1 = Region("NC_001139.9", 730000, 731500)
# reg2 = Region("NC_001139.9", 735000, 737000)    # INV   can't be detected by depth
# bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/aln.sort.bam"
# print(is_pair(reg1, reg2, bam_in))

def check_bed():
    import time
    time1 = time.time()
    bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/candidate_regions/candidate.bed"
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/aln.sort.bam"
    reg_ls  = []
    with open(bed_in, 'r') as f:
        for line in f:
            if line.strip().split()[3] != "clip_reg":
                continue
            reg_ls.append(Region(line.strip().split()[0], int(line.strip().split()[1]), int(line.strip().split()[2])))
    # bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in + ".bai")
    
    pre_reg = None
    cnt = 0
    # cnt2 = 0
    pairs_reg = []
    dup_reg = []
    pairs_merge_reg = []
    for reg in reg_ls:
        if not pre_reg:
            pre_reg = reg
        else:
            if is_pair(pre_reg, reg, bam_in):
                # print("pairs:{},{}-{} and {}-{}".format(pre_reg.chr_id, pre_reg.start, pre_reg.end, reg.start, reg.end))
                cnt += 1

                if pre_reg in pairs_reg:
                    dup_reg.append(pre_reg)
                pairs_reg.append(pre_reg)
                pairs_reg.append(reg)
                pairs_merge_reg.append(Region(pre_reg.chr_id, pre_reg.start, reg.end))
                
                
            pre_reg = reg
    print("total pairs:{}".format(cnt))
    for reg in pairs_merge_reg:
        print("{}:{}-{}".format(reg.chr_id, reg.start, reg.end))

    print("dup_reg_num:{}".format(len(dup_reg)))
    for reg in dup_reg:
        print("{}:{}-{}".format(reg.chr_id, reg.start, reg.end))
        
    time2 = time.time()
    print("cost time:{}s".format(time2 - time1))
# check_bed()



def check_flank(bedinfo_ls, i):
    '''intervals'''
    target_reg = [bedinfo_ls[i][0], bedinfo_ls[i][1], bedinfo_ls[i][2]]  # target low_dep reg
    merge_ls = []
    dis_bias = 100
    merge_ls.append(target_reg)
    ## 向前找
    for j in range(i-1, -1, -1):
        if target_reg[1] - bedinfo_ls[j][2] <= dis_bias:
            merge_ls.append([bedinfo_ls[j][0], bedinfo_ls[j][1], bedinfo_ls[j][2]])
        else:
            break
    ## 向后找
    for j in range(i+1, len(bedinfo_ls)):
        if bedinfo_ls[j][1] - target_reg[2] <= dis_bias:
            merge_ls.append([bedinfo_ls[j][0], bedinfo_ls[j][1], bedinfo_ls[j][2]])
        else:
            break
    return merge_ls

def clip_cov_reg_merge(bedinfo_ls_in):  # ctg bedibnfo_ls from candidate bed
    '''
    specify ctg, merge low_dep reg and flank clip_reg
    bed format: ctg, start, end, info(clip_reg、low_dep)
    遍历区间, 找到low_dep就往左右搜索, 看是否能够合并
    '''
    if len(bedinfo_ls_in) <= 1:
        return bedinfo_ls_in
    merge_ls = []   # 
    bedinfo_ls_out = []
    for i in range(len(bedinfo_ls_in)):
        if bedinfo_ls_in[i][3] == 'low_dep':
            res = check_flank(bedinfo_ls_in, i)
            if len(res) > 1:
                res.sort(key=lambda i:i[1])
                r_bound = max([i[2] for i in res])  # res[-1][2]
                merge_ls.extend(res)
                bedinfo_ls_out.append([res[0][0], res[0][1], r_bound, "low_dep_clip"])
                # print("merge regions from:", res, " To: {}:{}-{}".format(res[0][0], res[0][1], r_bound))
    ## 
    for bedinfo in bedinfo_ls_in:
        if [bedinfo[0],bedinfo[1],bedinfo[2]] in merge_ls:
            continue
        else:
            bedinfo_ls_out.append(bedinfo)
    # for 
    return bedinfo_ls_out

def check_merge():
    bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/candidate_regions/candidate.bed"
    pre_ctg = None
    ctg_bedinfo_ls = []
    ctg_bedinfo_merged_ls = []
    with open(bed_in,'r') as f:
        for line in f:
            filelds = line.strip().split('\t')
            if filelds[0] != pre_ctg:
                if len(ctg_bedinfo_ls) > 0:
                    print("\n--------------------{}--------------------".format(pre_ctg))
                    ctg_bedinfo_merged_ls = clip_cov_reg_merge(ctg_bedinfo_ls)
                    print(ctg_bedinfo_ls,"\n", ctg_bedinfo_merged_ls)
                    ctg_bedinfo_ls = []
                pre_ctg = filelds[0]
                ctg_bedinfo_ls.append([filelds[0], int(filelds[1]), int(filelds[2]), filelds[3]])
            else:
                ctg_bedinfo_ls.append([filelds[0], int(filelds[1]), int(filelds[2]), filelds[3]])
    
    # return ctg_bedinfo_merged_ls
# check_merge()



def is_pair2(reg1, reg2, bam_in):    # 一前一后两个区间，判断是否为一对pair。根据bam来判断
    if reg1.chr_id!= reg2.chr_id:
        return False
    max_connect_len = 20000
    if reg1.start - reg2.end > max_connect_len:
        return False
    ## 根据区间内的信息判断
    MIN_SUPPOPRT_NUM = 10
    MIN_MAPPING_QUALITY = 20
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in + ".bai")
    reads_set1 = set()
    reads_set2 = set()
    for read in bam_reader.fetch(reg1.chr_id, reg1.start, reg1.end):
        if read.is_secondary and read.mapping_quality < MIN_MAPPING_QUALITY:
            continue
        reads_set1.add(read.query_name)
    for read in bam_reader.fetch(reg1.chr_id, reg2.start, reg2.end):
        if read.is_secondary and read.mapping_quality < MIN_MAPPING_QUALITY:
            continue
        reads_set2.add(read.query_name)
    commons = reads_set1.intersection(reads_set2)
    print(commons)
    if len(commons) >= MIN_SUPPOPRT_NUM:
        return True
    return False

def clip_pair_cluster(reg_in, bam_in):   # cluster clip pair  ????
    bam_index = bam_in + ".bai"
    if len(reg_in) <= 1:
        return reg_in
    reg_out = []
    ctg = reg_in[0].chr_id
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index)
    max_connect_len = 20000
    reg_start = -1
    reg_end = -1
    for reg in reg_in:
        if reg.start - reg_end <= max_connect_len:  # 可能存在链接，能够配对
            pass
        else:
            if reg_start > -1:   # 非首次
                reg_out.append(Region(ctg, reg_start, reg_end))
                # reg_start = reg.start
                # reg_end = reg.end
            reg_start = reg.start
            reg_end = reg.end
    # for read in bam_reader.fetch():
    return reg_out


def test_samtools_cov():
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/aln.sort.bam"
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    reg = "NC_001144.5:470000-470200"
    # samtools coverage -m -l 1000 -q 20 -r NC_001144.5:460000-490000 aln.sort.bam
    res = pysam.coverage("-l", "1000", "-q", "20", "-r", reg, bam_in)
    print(res)
    for line in res.split("\n"):
        if not line or line.startswith("#"):
            continue
        print(line)
        
# test_samtools_cov()
# print(bool("a"))

def test():  
    """测试list和numpy数组的运算速度"""
    import time
    import numpy as np
    
    array_list = [1] * 100000000 
    print(len(array_list))
    
    t1 = time.time() 
    for i in range(0, len(array_list)-100, 100):
        cal_sum(array_list, i, i+100)  
    t2 = time.time()
    print("List elapsed time:%f"% (t2 - t1))
    
    t3 = time.time()
    array_np = np.array(array_list)
    t4 = time.time()
    print("Convert to Numpy time: %f" % (t4 - t3))
    
    t5 = time.time()
    for i in range(0, len(array_list)-100, 100):
        array_np[i:i+100].sum()
    t6 = time.time()
    print("Numpy elapsed time:%f"% (t6 - t5))

# test() 


#################################################
def check_flank(bedinfo_ls, i):
    '''intervals'''
    target_reg = [bedinfo_ls[i][0], bedinfo_ls[i][1], bedinfo_ls[i][2]]  # target low_dep reg
    merge_ls = []
    dis_bias = 100
    merge_ls.append(target_reg)
    ## 向前找
    for j in range(i-1, -1, -1):
        if target_reg[1] - bedinfo_ls[j][2] <= dis_bias:
            merge_ls.append([bedinfo_ls[j][0], bedinfo_ls[j][1], bedinfo_ls[j][2]])
        else:
            break
    ## 向后找
    for j in range(i+1, len(bedinfo_ls)):
        if bedinfo_ls[j][1] - target_reg[2] <= dis_bias:
            merge_ls.append([bedinfo_ls[j][0], bedinfo_ls[j][1], bedinfo_ls[j][2]])
        else:
            break
    return merge_ls

def clip_cov_reg_merge(bedinfo_ls_in):  # ctg bedibnfo_ls from candidate bed
    '''
    specify ctg, merge low_dep reg and flank clip_reg
    bed format: ctg, start, end, info(clip_reg、low_dep)
    遍历区间, 找到low_dep就往左右搜索, 看是否能够合并
    '''
    if len(bedinfo_ls_in) <= 1:
        return bedinfo_ls_in
    merge_ls = []   # 
    bedinfo_ls_out = []
    for i in range(len(bedinfo_ls_in)):
        if bedinfo_ls_in[i][3] == 'low_dep':
            res = check_flank(bedinfo_ls_in, i)
            if len(res) > 1:
                res.sort(key=lambda i:i[1])
                r_bound = max([i[2] for i in res])  # res[-1][2]
                merge_ls.extend(res)
                bedinfo_ls_out.append([res[0][0], res[0][1], r_bound, "low_dep_clip"])
                # print("merge regions from:", res, " To: {}:{}-{}".format(res[0][0], res[0][1], r_bound))
    ## 
    for bedinfo in bedinfo_ls_in:
        if [bedinfo[0],bedinfo[1],bedinfo[2]] in merge_ls:
            continue
        else:
            bedinfo_ls_out.append(bedinfo)
    # for 
    return bedinfo_ls_out

def can_span_by_reads(ctg, start, end, bam_in):
    min_support_reads = 3
    MIN_MAPPING_QUALITY = 30
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    support_reads_ls = []
    for read in bam_reader.fetch(ctg, start, end):
        if read.is_unmapped or read.is_secondary or read.mapping_quality < MIN_MAPPING_QUALITY:
            continue
        if read.reference_start < start and read.reference_end > end:
            support_reads_ls.append(read)
    if len(set(support_reads_ls)) != len(support_reads_ls):
        print("{}:{}-{}  dup reads".format(ctg, start, end))
    if len(support_reads_ls) >= min_support_reads:
        return True
    else:
        return False

def clip_reg_merge(bedinfo_ls_in, bam_in):
    if len(bedinfo_ls_in) <= 1:
        return bedinfo_ls_in
    bedinfo_ls_out = []
    ## collect clip reg
    clipinfo_ls = []
    for rec in bedinfo_ls_in:
        if rec[3] == 'clip_reg':
            clipinfo_ls.append(rec)
        else:   # "low_dep"
            bedinfo_ls_out.append([rec[0], rec[1], rec[2], "lowdep_asm"])     # low dep 标记为denovo asm
    
    ## 判断一下clipreg是否可以横跨
    reads_span_candidate_ls = []
    asm_candidate_ls = []
    for rec in clipinfo_ls:
        if can_span_by_reads(rec[0], rec[1], rec[2], bam_in):
            reads_span_candidate_ls.append(rec)
            # clipinfo_ls.remove(rec) 
        else:
            asm_candidate_ls.append(rec)
    
    ##################################################################################
    ## 对能横跨的作merge，主要是del
    pre_rec = None
    merge_flag = False
    for rec in reads_span_candidate_ls:
        if pre_rec:
            # if is_pair(pre_rec, rec, bam_in):
            if is_pair(Region(pre_rec[0], pre_rec[1], pre_rec[2]), Region(rec[0], rec[1], rec[2]), bam_in):
                bedinfo_ls_out.append([pre_rec[0], pre_rec[1], rec[2], "del_reads"])
                merge_flag = True
            else:
                if not merge_flag:
                    bedinfo_ls_out.append([pre_rec[0], pre_rec[1], pre_rec[2], "other_reads"])
                merge_flag = False
            pre_rec = rec
        else:
            pre_rec = rec
    if pre_rec:
        if not merge_flag:
            bedinfo_ls_out.append([pre_rec[0], pre_rec[1], pre_rec[2], "other_reads"])


    ## 对clip reg 进行聚类 merge
    start = -1
    end = -1
    cluster_dis = 20000
    for rec in asm_candidate_ls:
        if start > -1:
            if rec[1] - end < cluster_dis:
                end = rec[2]
            else:
                bedinfo_ls_out.append([rec[0], start, end, "clip_asm"])
                start = rec[1]
                end = rec[2]
        else:
            start = rec[1]
            end = rec[2]
    if start > -1:
        bedinfo_ls_out.append([rec[0], start, end, "clip_asm"])
    bedinfo_ls_out.sort(key=lambda x: x[1])
    # print(bedinfo_ls_out)
    return bedinfo_ls_out

def load_candidate_reg_bed(candidate_bed):
    bedinfo_dic = dict()    # ctg:ctg_bedinfo_ls
    pre_ctg = ""
    ctg_bedinfo_ls = []
    with open(candidate_bed, 'r') as f:
        for line in f:
            fields = line.strip().split("\t")
            if fields[0] != pre_ctg:
                if pre_ctg:
                    bedinfo_dic[pre_ctg] = ctg_bedinfo_ls   # 
                    ctg_bedinfo_ls = []
                pre_ctg = fields[0]
                ctg_bedinfo_ls.append([fields[0], int(fields[1]), int(fields[2]), fields[3]])
            else:   
                ctg_bedinfo_ls.append([fields[0], int(fields[1]), int(fields[2]), fields[3]])
    if not ctg_bedinfo_ls:
        bedinfo_dic[pre_ctg] = ctg_bedinfo_ls
    return bedinfo_dic

def test():
    bed_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/candidate_regions2/candidate.bed"
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/aln.sort.bam"
    bedinfo_dic = load_candidate_reg_bed(bed_in)
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    # print(bedinfo_dic)
    read_cnt, lowdep_cnt, clip_cnt = 0, 0, 0
    asm_ls = []
    all_ls = []
    for ctg in bam_reader.references:
        print("\nProcess {}".format(ctg))
        if ctg not in bedinfo_dic.keys():
            continue
        ctg_bedinfo_ls = bedinfo_dic[ctg]
        ls = clip_cov_reg_merge(ctg_bedinfo_ls)
        ls = clip_reg_merge(ls, bam_in)
        all_ls.extend(ls)
        for rec in ls:
            print(rec)
            if rec[3].endswith("reads"):
                read_cnt += 1
            elif rec[3].endswith("asm"):
                if rec[3].startswith("lowdep"):
                    lowdep_cnt += 1
                else:
                    clip_cnt += 1
                asm_ls.append(Region(rec[0], rec[1], rec[2]))
    print(read_cnt, lowdep_cnt, clip_cnt)
    print("\n asm:")
    for reg in asm_ls:
        print(reg)
    bed_out = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_op/test.bed"
    fo = open(bed_out, "w")
    for rec in all_ls:
        fo.write("{}\t{}\t{}\t{}\n".format(rec[0], rec[1], rec[2], rec[3]))
    fo.close
    return
# test()
#################################################

# python /public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/create_consensus_by_bed/get_fasta_consensus2.py --fasta /public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/REF/yeast/GCF_000146045.2_R64_genomic.fna --bam /public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/aln.sort.bam --contig_ls NC_001133.9,NC_001134.8,NC_001135.5,NC_001136.10,NC_001137.3,NC_001138.5,NC_001139.9,NC_001140.6,NC_001141.2,NC_001142.9,NC_001143.9,NC_001144.5,NC_001145.3,NC_001146.8,NC_001147.6,NC_001148.4,NC_001224.1 --out-dir /public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/CRR198362/mypipe/test_pipe3 --bed candidate_regions2/candidate.bed --min-contig 0 -t 24


def cal_avg_depth2(bam_in, reg):  # cal depth of a contig, return list of depth, length == contig length
    # reg = reg.chr_id + ":" + str(reg.start) + "-" + str(reg.end)
    '''time samtools depth -a -Q $min_MQ -r $REG ${bam_file} > $out_prefix.depth    
    # 不加-J del处depth为0'''
    min_MQ = 40
    min_len = 5000
    # reg_len = int(reg.split(":")[2]) - int(reg.split(":")[1])
    depth_ls = []
    depth_stream = pysam.depth("-a", "-Q", str(min_MQ), "-l", str(min_len), "-g", "0x800", "-r", reg, bam_in)    # depth_stream接收stdout
    for line in depth_stream.split("\n"):
        line_ls = line.split()
        if not line_ls:
            continue
        depth_ls.append(int(line_ls[2]))
    return sum(depth_ls) / len(depth_ls)
# reg1 = "chrXI:75430-79067"
reg2 = "chrV:425500-429515"
print(cal_avg_depth2("/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/yeast/test_snakemake/snakemake_res/CRR198362_40X/step1_mapping/aln.sort.bam", reg2))

###########################################################################################################################################################################################################
def convert_reference_pos_to_raw_pos2(read, candidate_pos):  
    # 输入read信息和要求的参考上的位置，将参考的位置转换为read_querysequence上的位置
    '''注意:这个位置并不一定是读数原来的位置'''
    candidate_pos = set(candidate_pos)
    raw_ref_pos_map={}
    ref_pos = read.reference_start
    query_pos = 0
    for (ct,cl) in read.cigartuples:
        if ct==0:
            for i in range(cl):
                query_pos+=1
                ref_pos+=1
                if ref_pos in candidate_pos:
                    raw_ref_pos_map[ref_pos] = query_pos
        elif ct==1:
            query_pos+=cl
        elif ct==2:
            for i in range(cl):
                ref_pos+=1
                if ref_pos in candidate_pos:
                    raw_ref_pos_map[ref_pos] = query_pos
        elif ct==4:
            query_pos+=cl
        else:
            continue
    # return raw_ref_pos_map,read_length
    return raw_ref_pos_map

def solve_by_reads(ctg, start, end, support_ls):
    '''选取最优的一条读数序列，我们这里根据比对长度选择比对最长的一条序列'''
    if end <= start:
        print("{}:{}-{} ERROR REGION !!!".format(ctg, start, end))
        return ""
    solve_seq = ""
    ## 选择最佳读数
    max_len = -1
    best_read = None
    for read in support_ls:
        if read.query_alignment_length > max_len:
            max_len = read.query_alignment_length
            best_read = read
    ## 根据最佳读数选择   //读数reverse  不用作反转
    ref_to_read = convert_reference_pos_to_raw_pos2(best_read, [start, end])
    query_start = ref_to_read.get(start, -1)
    query_end = ref_to_read.get(end, -1)
    if query_start < 0 or query_end < 0 or query_start > query_end :
        print("{}:{}-{} ERROR!!!".format(ctg, start, end))
        return ""
    solve_seq = best_read.query_sequence[query_start:query_end]
    print("best_read:", best_read.query_name)
    return solve_seq

def judge_span_by_reads(ctg, start, end, bam_in):
    min_support_reads = 3
    MIN_MAPPING_QUALITY = 60
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in+".bai")
    support_reads_ls = []
    for read in bam_reader.fetch(ctg, start, end):
        if read.is_unmapped or read.is_secondary or read.mapping_quality < MIN_MAPPING_QUALITY:
            continue
        if read.reference_start < start-1000 and read.reference_end > end + 1000:
            support_reads_ls.append(read)
    if len(set(support_reads_ls)) != len(support_reads_ls):
        print("{}:{}-{}  dup reads".format(ctg, start, end))
    if len(support_reads_ls) >= min_support_reads:
        return solve_by_reads(ctg, start, end, support_reads_ls)
        # return solve
    return ""
def test3():
    bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/yeast/test_snakemake/snakemake_res/CRR198362_40X/step1_mapping/aln.sort.bam"
    ctg = "chrV"
    start = 433001
    end = 436000
    print(judge_span_by_reads(ctg, start, end, bam_in))
# test3()