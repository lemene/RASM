import yaml
import os
import math
import numpy as np
import gzip
import multiprocessing as mp
import subprocess
import pysam
import re
import time

import itertools

from functools import partial
from collections import namedtuple,defaultdict, Counter


from depth import *
from pileup import *
from clip import *

import utils


def cal_idels(read:pysam.AlignedRead, reg_start, reg_end):
    '''计算ref_start-ref_end之间的indels'''
    ref_pos = read.reference_start

    ins_num = 0
    del_num = 0
    sum_ins_num = 0
    sum_del_num = 0
    query_pos = 0
    match_num = 0
    min_size = 30
    
    rpos = read.reference_start
    qpos = 0
    for op, op_len in read.cigartuples:
        if ref_pos > reg_end: break
        if op == 0:     # match
            if op_len + ref_pos < reg_start or ref_pos > reg_end:
                pass
            elif op_len + ref_pos < reg_end:
                match_num += op_len
            else:
                match_num += reg_end - ref_pos
            query_pos += op_len
            ref_pos += op_len
        elif op == 1:   # ins
            query_pos += 1
            if ref_pos > reg_start and ref_pos < reg_end:
                # print("Add ins")
                if op_len > min_size:
                    ins_num += op_len
            sum_ins_num += op_len
        elif op == 2:   # del
            if op_len + ref_pos < reg_start:
                pass
            elif op_len + ref_pos < reg_end:    # []
                if op_len > min_size:
                    del_num += op_len
            else:   # > reg_end
                if op_len > min_size:
                    del_num += op_len
            ref_pos += op_len
            sum_del_num += op_len
        elif op == 4:   # soft clip
            query_pos += op_len
        else:
            continue
    # print(ins_num, del_num)
    return ins_num, del_num


class Detector:
    def __init__(self, cfg, wrkdir, bam_fname, asm_fname):
        self.cfg = yaml.safe_load(open(cfg))
        self.wrkdir = wrkdir
        self.bam_fname = bam_fname
        self.asm_fname = asm_fname
        self.infos = {}
        
        utils.safe_make_dir(wrkdir)
        utils.enable_logging(os.path.join(self.wrkdir, "lrmd.log"))
        utils.logger.info("Start detecting misassemblies")

        utils.logger.info("Config参数设置: {}".format(self.cfg))
        
        
    def get_contigs(self, min_contig):
        bam = pysam.AlignmentFile(self.bam_fname, "rb")
        ctgs = [(ref, rlen) for ref, rlen in zip(bam.references, bam.lengths) if rlen >= min_contig]
        return ctgs
    

    
    def detect(self, threads, min_contig):

        bam = pysam.AlignmentFile(self.bam_fname, "rb")
        ctgs = [(ref, rlen) for ref, rlen in zip(bam.references, bam.lengths) if rlen >= min_contig and ref == "CP132235.1"]

        candidates = []

        win_size = self.cfg["win_size"]
        stride = self.cfg["stride"]
        depth_min_ratio = self.cfg["depth_min_ratio"]
        depth_max_ratio = self.cfg["depth_max_ratio"]
        min_clip_len = self.cfg["min_clip_len"]
        min_mapq = self.cfg["min_mapq"]
        min_clip_num = self.cfg["min_clip_num"]
        max_diff_ratio = self.cfg["max_diff_ratio"]
        
        if self.cfg['apply_depth']:
            depth_info = DepthInfo(self.wrkdir)
            depth_info.build(ctgs, self.bam_fname, threads)
            candidates.extend(depth_info.get_mis_candidates( win_size, stride, depth_min_ratio, depth_max_ratio))
            utils.logger.info("+ depth feature: %d" % len(candidates))

            
        if self.cfg["apply_clip"]:
            clip_info = ClipInfo()
            clip_info.build(ctgs, self.bam_fname, threads, min_mapq, min_clip_len)
            candidates.extend(clip_info.get_mis_candidates(win_size, stride, min_clip_num))
            utils.logger.info("+ clip feature: %d" % len(candidates))

        if self.cfg["apply_pileup"]:
            pileup_info = PileupInfo(self.wrkdir)
            pileup_info.build(ctgs, self.bam_fname, self.asm_fname, threads, min_mapq)
            candidates.extend(pileup_info.get_mis_candidates(win_size, stride, max_diff_ratio))
            utils.logger.info("+ pileup feature: %d" % len(candidates))


        merged = self.merge_segments(candidates, 5000)
        
        for ctg_name, start, end in merged:
            utils.logger.info("merged: %s:%d-%d" % (ctg_name, start, end))
        utils.logger.info("merged candidate size = %d", len(merged))
        
        self.filter(merged)

    def run_find_pipe(self, ref, bam, ctgs, out_dir, threads):
        self.depth_infos.build(ctgs, bam, threads)

        ## ----------------1、find candidate----------------
        # dirs
        candidate_dir = os.path.join(out_dir, "candidate")
        if not os.path.isdir(candidate_dir):os.mkdir(candidate_dir)
        info_bed_dir = os.path.join(out_dir, "info")
        if not os.path.isdir(info_bed_dir):os.mkdir(info_bed_dir)
        filtered2_dir = os.path.join(out_dir, "filtered2")
        if not os.path.isdir(filtered2_dir):os.mkdir(filtered2_dir)
        pileup_dir = os.path.join(out_dir, "pileup")
        if not os.path.isdir(pileup_dir):os.mkdir(pileup_dir)

        # run find
        pool = mp.Pool(processes=threads)
        for ctg in ctgs:
            pool.apply_async(partial(self.find_candidate), args=(ctg[0], ctg[1], ref, bam, out_dir, self.cfg["step2"]))
        
        pool.close() 
        pool.join()

    def find_candidate(self, ctg, reg_end, ref, bam, out_dir, config):
        bam_reader = pysam.AlignmentFile(bam, "rb")
        
        reg_start = 0
        ctg_len = reg_end
        # OUT: 
        # 
        reg_len = reg_end - reg_start
        win_size = config["win_size"]
        stride = win_size // 2
        win_num = (reg_len - win_size + stride - 1) // stride + 1 
        
        # clip params
        clip_params = config["clip_params"]
        min_clip_portion = clip_params["min_clip_portion"]
        min_clip_num = min_clip_portion * self.depth_infos.median
        min_clip_len = clip_params["min_clip_len"]

        # pileup params
        pileup_params = config["pileup_params"]
        win_size, step_size, min_correct_portion, max_differ_portion, max_disagree_portion, cluster_dis \
            = pileup_params["win_size"], pileup_params["step_size"], pileup_params["min_correct_portion"], \
            pileup_params["max_differ_portion"], pileup_params["max_disagree_portion"], pileup_params["cluster_dis"]
        
        ## -----------------------一、cal info-----------------------

        win_clip_num = self.get_clip_info(ctg, clip_params, bam_reader, win_num, stride, win_size, min_clip_len)
        
        candidates = []

        if config['apply_dp']:
            candidates.extend(self.depth_infos.get_mis_candidates(ctg, win_size, stride, config["depth_min_ratio"], config["depth_max_ratio"]))

        if config['apply_pileup']:
            pipeup_info = PileupInfo(self.wrkdir)
            pipeup_info.build(bam, ctg, ref, 10)
            candidates.extend(pipeup_info.get_mis_candidates(   ))

        if config['apply_clip']:
            cilp_info = ClipInfo()
            clip_info.build()


        for i in range(win_num):
            s, e = i * stride, min(i * stride + win_size, ctg_len)

            if config['apply_clip'] and (win_clip_num[i] >= min_clip_num) or \
               config['apply_pileup'] and (win_correct_portion[i] <= min_correct_portion or win_differ_portion[i] >= max_differ_portion or win_disagree_portion[i] >= max_disagree_portion):
                
                candidate_ls.append([ctg, s, e])

        merged = self.merge_segments(candidates, 5000)
        
        self.filter2(bam, ctg, ctg_len, merged, out_dir, config)

    def is_read_trivial(self, read, mapq):
        return read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < mapq

    def is_covering_region(self, ctg, start, end, bam, ref):
        check_win_size = self.cfg["check_win_size"]
        utils.logger.info("filter: %s:%d-%d" % (ctg, start, end))
        win_iter = utils.WinIterator(end - start, check_win_size, check_win_size)
        for s, e in win_iter:
            s += start
            e += start

            if not self.is_coverage_win(ctg, s, e, bam, ref):
                return False
        return True

    def calc_distance(self, start, end, read, ref):
        rpos = 0 # read.reference_start
        rseq = ref.fetch(read.reference_name, read.reference_start, read.reference_end)
        distance = np.zeros(len(rseq))
        qpos = 0
        qseq = read.query_sequence
        for op, op_len in read.cigartuples:
            #utils.logger.info("op: %d %d %d %d" % (op, op_len, rpos+read.reference_start, qpos))
            if op == 0:         # M
                #for i in range(op_len):
                #    utils.logger.info("%s %s\n", rseq[rpos+i], qseq[qpos+i])
                qpos += op_len
                rpos += op_len
                pass
            elif op == 1:   # I insertion
                distance[rpos] += op_len
                qpos += op_len
            elif op == 2:   # D deletion
                distance[rpos:rpos+op_len] += 1
                rpos += op_len
            elif op == 7:   # =
                qpos += op_len
                rpos += op_len
            elif op == 8:   # X
                distance[rpos:rpos+op_len] += 1
                qpos += op_len
                rpos += op_len
            else:
                pass        # 不处理
        
        #utils.logger.info("distance: %s" % (list(distance)))
        #utils.logger.info("distance: %d %d %d %d" % (sum(distance), sum(distance[start-read.reference_start:end-read.reference_start]), start, end))
        #assert 0
       
        return sum(distance[start-read.reference_start:end-read.reference_start]) / (end-start)
    

    def stat_cigar(self, start, end, read, ref):
        
        pass

    def is_coverage_win(self, ctg, start, end, bam, ref):
        
        min_mapq = self.cfg["min_mapq"]
        min_clip_len = self.cfg["min_clip_len"]
        min_distance = self.cfg["min_distance"]

        span = []
        ctg_len = bam.get_reference_length(ctg)
        for read in bam.fetch(ctg, start, end):
            utils.logger.info("ccc %s %d %d %d %d" % (read.qname, read.mapping_quality, read.is_unmapped, read.is_secondary, read.is_supplementary))
            if (self.is_read_trivial(read, min_mapq)):
                continue

            utils.logger.info("cigar %s %d %d %d %d" % (read.qname, read.cigartuples[0][0], read.cigartuples[0][1], 
                                                        read.cigartuples[-1][0], read.cigartuples[-1][1]))
            cigar = read.cigartuples
            left = read.cigartuples[0]
            if left[0] == 4 or left[0] == 5:
                if min(left[1], read.reference_start) > min_clip_len:
                    continue
            right = read.cigartuples[-1]
            if right[0] == 4 or right[0] == 5:
                if min(right[1], ctg_len - read.reference_end) > min_clip_len:
                    continue
            
            utils.logger.info("distance: %s %f, %d", read.qname, self.calc_distance(start, end, read, ref), min_clip_len)
            if self.calc_distance(start, end, read, ref) < 0.2:
                
                if read.reference_start <= max(0, start - 500) and read.reference_end >= min(ctg_len, end + 500):
                    span.append(read)
                    utils.logger.info("span: %d", len(span))

        
        return len(span) >= 3
    

    def win_check(self, ctg, start, end, ctg_len, bam_reader, params):
        print("Check:{}:{}-{}".format(ctg, start, end))
        '''
        True: mis
        False: no mis
        1、对深度的限制，一些重复区不一定能从cigar检测出mis。或者存在干扰
        2、对
        hifi使用不同的参数
        '''
        min_span_num = params['min_span_num']
        min_supp_portion = params['min_supp_portion']  # 有多少成的读数表明没有mis
        min_MQ = params['min_MQ']    # 10    map_Q太严了，导致大量的比对过滤掉了->3
        MIN_ALIGN_LENGTH = params['MIN_ALIGN_LENGTH']
        MIN_ALIGN_RATE = params['MIN_ALIGN_RATE']
        ins_threshold = params['ins_threshold']
        del_threshold = params['del_threshold']
        indel_threshold = ins_threshold + del_threshold
        min_clip_len = params['min_clip_len']
        bound_span = params['bound_span']   # 5000/10000
        # 

        span_ls = []
        for read in bam_reader.fetch(ctg, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < min_MQ:
                continue
            # for bound
            is_bound = False
            if read.reference_start <= max(500, bound_span) or read.reference_end >= min(ctg_len - 500, ctg_len - bound_span):
                MIN_ALIGN_LENGTH = MIN_ALIGN_LENGTH - 2000
                MIN_ALIGN_RATE = MIN_ALIGN_RATE - 0.1
                is_bound = True
            ## 
            if read.query_alignment_length < MIN_ALIGN_LENGTH or (read.query_alignment_length / read.query_length) < MIN_ALIGN_RATE:
                continue
            cigar = read.cigartuples
            left = cigar[0]
            if left[0] == 4 or left[0] == 5:
                if left[1] > min_clip_len and not is_bound:    # left clip非左边界，且其长度大于阈值
                    continue
            right = cigar[-1]
            if right[0] == 4 or right[0] == 5:
                if right[1] > min_clip_len and not is_bound: # right clip非右边界，且其长度大于阈值
                    continue
            if read.reference_start <= max(0, start - 500) and read.reference_end >= min(ctg_len, end + 500):
                span_ls.append(read)
        
        #  
        if len(span_ls) < min_span_num and start >= 5000 and end <= ctg_len - 5000:
            return True, "low_supp_mis"
        else:
            min_span_num = max(min_span_num, math.ceil(min_supp_portion*len(span_ls)))
            # cal indel of the span
            Inss_ls = []
            Dels_ls = []
            Indels_ls = []
            for read in span_ls:   # 收集span_ls读数的indels
                read_l = max(start - 5000, read.reference_start + 100)
                read_r = min(end + 5000, read.reference_end - 100)
                indels = cal_idels(read, read_l, read_r)
                Indels_ls.append(indels[0] + indels[1])
                Inss_ls.append(indels[0])
                Dels_ls.append(indels[1])
            # 
            Inss_ls.sort()
            Dels_ls.sort()
            Indels_ls.sort()
            # check indels
            avg_ins = sum(Inss_ls[:min_span_num]) // min_span_num
            avg_del = sum(Dels_ls[:min_span_num]) // min_span_num
            avg_indels = sum(Indels_ls[:min_span_num]) // min_span_num
            # print(Inss_ls, Dels_ls, Indels_ls)
            if avg_ins < ins_threshold and avg_del < del_threshold and avg_indels < indel_threshold: # no mis
                return False, "No_mis"
            else:
                return True, "Reads_mis"

    def filter(self, candidates):
        filtered = []

        bam = pysam.AlignmentFile(self.bam_fname, "rb")
        ref = pysam.FastaFile(self.asm_fname)
        for ctg_name, start, end in candidates:
            #if start != 3114800: continue
            if not self.is_covering_region(ctg_name, start, end, bam, ref):
                filtered.append((ctg_name, start, end))
        
        utils.logger.info("misassembly size = {}".format(len(filtered)))
        
        
        with open(os.path.join(self.wrkdir, "misassembly.bed"), "w") as f:
            for ctg_name, start, end in filtered:
                f.write("%s\t%d\t%d\n" % (ctg_name, start, end))


    def get_pileup_info(self, pileup_params, ctg, reg_start, reg_end, out_dir, ref, bam, win_size, step_size):    
        min_MQ = pileup_params['min_MQ']

        region = ctg + ":" + str(0) + "-" + str(reg_end)
        pileup_dir = os.path.join(out_dir, "pileup")
        pileup_file = os.path.join(pileup_dir, region + ".pileup.txt")
        pileup_cmd = ["samtools mpileup -B", "-q", str(min_MQ), "-aa", "-d 200","-r", region, "-f", ref, bam, "-o", pileup_file]
        subprocess.check_call(" ".join(pileup_cmd), shell=True)

        pileup_dict={"contig":[],"correct":[],"ambiguous":[],"insert":[],"deletion":[],
                    "disagree":[],"depth":[],"differ":[]}
        # for line in pileup_stream.split("\n"):
        with open(pileup_file, "r") as f:
            for line in f:
                record = line.split()
                assert len(record) == 6

                ## 
                pileup_dict['contig'].append(record[0])
                match_detail=record[4]
                pileup_dict['correct'].append(match_detail.count('.')+match_detail.count(','))
                pileup_dict["depth"].append(int(record[3]))     # 反映的是pileup计算的depth，与实际计算的可以不一致
                st = ''.join(re.split('[\+|\-][0-9]+[ATCGatcg]+', match_detail))
                st_counter = Counter(st)
                disagree_numd = st_counter['a'] + st_counter['A'] + st_counter['g'] + st_counter['G'] + st_counter['c'] + st_counter['C'] + st_counter['t'] + st_counter['T']   # 单碱基的不一致，snp
                pileup_dict["disagree"].append(disagree_numd)
                match_counter = Counter(match_detail)
                differ_numd = match_counter['a'] + match_counter['A'] + match_counter['g'] + match_counter['G'] + match_counter['c'] + match_counter['C'] + match_counter['t'] + match_counter['T'] - disagree_numd # 包括插入删除的大小信息
                # ins_numd = 
                pileup_dict["differ"].append(differ_numd)
        
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
                continue
            window_pileup_dict["correct_portion"].append(np.sum(pileup_dict['correct'][start:end])/total)
            window_pileup_dict["disagree_portion"].append(np.sum(pileup_dict["disagree"][start:end])/total)
            window_pileup_dict["differ_portion"].append(np.sum(pileup_dict["differ"][start:end])/total) 

        win_correct_portion = window_pileup_dict["correct_portion"]
        win_differ_portion = window_pileup_dict["differ_portion"]
        win_disagree_portion = window_pileup_dict["disagree_portion"]
        return win_correct_portion,  win_differ_portion, win_disagree_portion


    def get_clip_info(self, ctg, clip_params, bam_reader, win_num, stride, win_size, min_clip_len):

        ctg_len = bam_reader.get_reference_length(ctg)
        clip_ls = [0] * ctg_len 

        min_MQ = clip_params['min_MQ']
        for read in bam_reader.fetch(ctg):
            if read.is_unmapped or read.is_secondary or read.mapping_quality < min_MQ:
                continue

            left = read.cigartuples[0]
            if (left[0] == 4 or left[0] == 5)  and left[1] >= min_clip_len:
                    clip_ls[read.reference_start] += 1

            right = read.cigartuples[-1]
            if (right[0] == 4 or right[0] == 5) and right[1] >= min_clip_len:
                    clip_ls[read.reference_end-1] += 1
                    
        win_clip_num = [0] * win_num
        for i in range(win_num):   
            win_clip_num[i] = sum(clip_ls[i*stride: min(i*stride+win_size, ctg_len)])
        
        return win_clip_num
    
    
    def merge_segments(self, segs, distance):
        '''segs: [(ctg, start, end)]'''
        
        assert len(segs) >= 0
        segs.sort()             # 字典順序
        merged = [list(segs[0])]

        for s in segs[1:]:
            if s[0] == merged[-1][0] and s[1] - merged[-1][2] <= distance:
                merged[-1][2] = max(merged[-1][2], s[2])
            else:
                merged.append(list(s))
        return merged


def cal_Ne(work_dir, ctg, ctg_len, min_dp, mm_rate):
    Ne, Na = 0, 0
    
    region = ctg + ":" + str(0) + "-" + str(ctg_len)
    pileup_file = os.path.join(work_dir, "pileup", region + ".pileup.txt")
    mis_bed = os.path.join(work_dir, "filtered2", region + ".bed")
    mis_ls = []
    with open(mis_bed, "r") as f:
        for line in f:
            fields = line.strip().split()
            mis_ls.append([fields[0], int(fields[1]), int(fields[2])])        
    if len(mis_ls) > 0:
        mis = mis_ls[0]
    else:            
        mis = [ctg, -1, -1]
    mis_length = 0
    mis_length += mis[2] - mis[1]
    mis_idx = 0
    
    with open(pileup_file, "r") as f:
        for line in f:
            items = line.strip().split("\t")
            pos = int(items[1])
            if pos >= mis[1] and pos <= mis[2]: 
                continue 

            if pos > mis[2]:  # over this mis
                if mis_idx + 1 < len(mis_ls):
                    mis_idx += 1
                    mis = mis_ls[mis_idx]
                    mis_length += mis[2] - mis[1]

            depth = int(items[3])
            match = items[4].count('.') + items[4].count(',')
            if depth >= min_dp:
                Na += 1
                if match / depth < mm_rate:
                    Ne += 1
    return Ne, Na

def get_qv(work_dir, ctgs, threads):

    simple_file = os.path.join(work_dir, "simple.ststs")
    f = open(simple_file, "a+")
    min_dp = 5
    mm_rate = 0.7   # 根据读数
    # 
    results = []
    pool = mp.Pool(processes=threads)
    for ctg, ctglen in ctgs:
        results.append(pool.apply_async(cal_Ne, args=(work_dir, ctg, ctglen, min_dp, mm_rate)))

    pool.close() 
    pool.join()

    Ne = 0
    Na = 0
    for res in results:
        ctg_Ne, ctg_Na = res.get()
        Ne += ctg_Ne
        Na += ctg_Na
    qv = -10 * math.log10(Ne/Na)

    f.write("QV: {}\n".format(qv))
    f.close()