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
from summary import *

import utils


class Detector:
    def __init__(self, cfg, wrkdir, bam_fname, asm_fname):
        self.cfg = yaml.safe_load(open(cfg))
        self.wrkdir = wrkdir
        self.bam_fname = bam_fname
        self.asm_fname = asm_fname
        self.infos = {}
        
        utils.safe_make_dir(wrkdir)
        utils.logger.info("Start detecting misassemblies")

        utils.logger.info("Config参数设置: {}".format(self.cfg))

        self.summary = Summary(bam_fname, asm_fname)
        
        
    def get_contigs(self, min_contig):
        bam = pysam.AlignmentFile(self.bam_fname, "rb")
        ctgs = [(ref, rlen) for ref, rlen in zip(bam.references, bam.lengths) if rlen >= min_contig]
        return ctgs
    

    
    def detect(self, threads, min_contig):
        bam = pysam.AlignmentFile(self.bam_fname, "rb")
        ctgs = [(ref, rlen) for ref, rlen in zip(bam.references, bam.lengths) if rlen >= min_contig]

        candidates = []

        win_size = self.cfg["win_size"]
        stride = self.cfg["stride"]
        min_clip_len = self.cfg["min_clip_len"]
        min_mapq = self.cfg["min_mapq"]
        min_clip_num = self.cfg["min_clip_num"]
        max_diff_ratio = self.cfg["max_diff_ratio"]
        
        if self.cfg['apply_depth']:
            depth_info = DepthInfo()
            depth_info.build(ctgs, self.summary, threads)
            cands = depth_info.get_mis_candidates(win_size, stride)
            candidates.extend(cands)
            utils.logger.info("+ depth feature: %d" % len(candidates))

        if self.cfg["apply_clip"]:
            clip_info = ClipInfo()
            clip_info.build(ctgs, self.bam_fname, threads, min_clip_len)
            candidates.extend(clip_info.get_mis_candidates(win_size, stride, min_clip_num))
            utils.logger.info("+ clip feature: %d" % len(candidates))

        if self.cfg["apply_pileup"]:
            pileup_info = PileupInfo()
            pileup_info.build(ctgs, self.summary, threads)
            candidates.extend(pileup_info.get_mis_candidates(win_size, stride, max_diff_ratio))
            utils.logger.info("+ pileup feature: %d" % len(candidates))


        
        merged = self.merge_segments(candidates, 5000) if len(candidates) > 0 else []
        
        for ctg_name, start, end in merged:
            utils.logger.info("merged: %s:%d-%d" % (ctg_name, start, end))
        utils.logger.info("merged candidate size = %d", len(merged))
        
        misassembly = merged
        misassembly = self.filter(merged)

        with open(os.path.join(self.wrkdir, "misassembly.bed"), "w") as f:
            for ctg_name, start, end in misassembly:
                f.write("%s\t%d\t%d\n" % (ctg_name, start, end))

    def is_read_trivial(self, read, mapq):
        return read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < mapq

    def is_covering_region(self, ctg, start, end, bam, ref):
        check_win_size = self.cfg["check_win_size"]
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
            if op == 0:         # M
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
        
       
        return sum(distance[start-read.reference_start:end-read.reference_start]) / (end-start)
    


    def is_coverage_win(self, ctg, start, end, bam, ref):
        
        min_mapq = self.cfg["min_mapq"]
        min_clip_len = self.cfg["min_clip_len"]
        min_distance = self.cfg["min_distance"]

        span = []
        clip_num = 0
        lowqual_num = 0
        ctg_len = bam.get_reference_length(ctg)
        for read in bam.fetch(ctg, start, end):
            if (self.is_read_trivial(read, min_mapq)):
                continue

            cigar = read.cigartuples
            left = read.cigartuples[0]
            if left[0] == 4 or left[0] == 5:
                if min(left[1], read.reference_start) > min_clip_len:
                    if read.reference_start >= start and  read.reference_start < end:
                        clip_num += 1
                    continue
            right = read.cigartuples[-1]
            if right[0] == 4 or right[0] == 5:
                if min(right[1], ctg_len - read.reference_end) > min_clip_len:
                    if read.reference_end > start and  read.reference_end <= end:
                        clip_num += 1
                    continue
            
            if self.calc_distance(start, end, read, ref) < 0.2:
                if read.reference_start <= max(0, start - 500) and read.reference_end >= min(ctg_len, end + 500):
                    span.append(read)
            else:
                lowqual_num += 1
        
        print("is_coverage_win: %s:%d-%d %d %d %d" % (ctg, start, end, len(span), lowqual_num, clip_num))
        return len(span) >= 3 and lowqual_num < min(10, len(span)/2) and clip_num < min(10, len(span)/2) 
        ##return len(span) >= 3
    

    def filter(self, candidates):
        filtered = []

        bam = pysam.AlignmentFile(self.bam_fname, "rb")
        ref = pysam.FastaFile(self.asm_fname)
        for ctg_name, start, end in candidates:
            if not self.is_covering_region(ctg_name, start, end, bam, ref):
                filtered.append((ctg_name, start, end))
        
        utils.logger.info("misassembly size = {}".format(len(filtered)))
        return filtered
        

 
    def merge_segments(self, segs, distance):
        '''segs: [(ctg, start, end)]'''
        
        assert len(segs) > 0
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