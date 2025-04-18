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
from candidate import *

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

        self.cfg["threads"] = threads
        
    
        candidates = Candidate().get_mis_candidates(ctgs, self.summary, self.cfg)
        utils.logger.info(f"Find candidate {len(candidates)}")
        merged = self.merge_segments(candidates, 5000) if len(candidates) > 0 else []
        
        for ctg_name, start, end in merged:
            utils.logger.debug("merged: %s:%d-%d" % (ctg_name, start, end))
        utils.logger.info("merged candidate size = %d", len(merged))
        
        misassemblies = self.verify_candidates(merged)

        self.dump_misassemblies(os.path.join(self.wrkdir, "misassembly.bed"), misassemblies)


    def is_read_trivial(self, read, mapq):
        return read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < mapq

    def is_covering_region(self, ctg, start, end, bam, ref):
        check_win_size = self.cfg["check_win_size"]
        (lcov, rcov) = self.check_surrounding_region(ctg, start, end)
        th_cov = max(3, min(lcov, rcov) / 2)
        win_iter = utils.WinIterator(end - start, check_win_size, check_win_size)
        for s, e in win_iter:
            s += start
            e += start

            if not self.is_coverage_win(ctg, s, e, bam, ref, th_cov):
                return False
        return True
    
    def check_surrounding_region(self, ctg, start, end):
        min_mapq = self.cfg["min_mapq"]
        min_clip_len = self.cfg["min_clip_len"]
        min_distance = self.cfg["min_distance"]

        start0 = max(0, start - 1000)
        end1 = min(end+1000, self.summary.get_contig_length(ctg))

        # 计算平均覆盖度
        lcov = -1   # 表示左边没有区域
        if start0 + 100 < start:    # 认为超过100才有统计意义
            lcov = self.summary.get_average_coverage(ctg, start0, start)
        rcov = -1
        if end + 100 < end1:
            rcov = self.summary.get_average_coverage(ctg, end, end1)
        return lcov, rcov
        


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
    


    def is_coverage_win(self, ctg, start, end, bam, ref, th_cov):
        
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
        
        print(f"is_coverage_win: {ctg}:{start}-{end} {len(span)} >= {th_cov} {lowqual_num} {clip_num}")
        return len(span) >= th_cov
    

    def verify_candidates(self, candidates):
        verified = []

        bam = pysam.AlignmentFile(self.bam_fname, "rb")
        ref = pysam.FastaFile(self.asm_fname)
        for ctg_name, start, end in candidates:
            if not self.is_covering_region(ctg_name, start, end, bam, ref):
                verified.append((ctg_name, start, end))
        
        utils.logger.info("misassembly size = {}".format(len(verified)))
        return verified
        

 
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

    def dump_misassemblies(self, fname, misassembies):
        with open(os.path.join(self.wrkdir, "misassembly.bed"), "w") as f:
            for ctg_name, start, end in misassembies:
                f.write("%s\t%d\t%d\n" % (ctg_name, start, end))
