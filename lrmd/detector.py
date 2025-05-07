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


from summary import *
from candidate import *

import utils
from enum import Enum

class RegionType(Enum):
    NORMAL = 0      # 正常区域
    ERROR = 1       # 组装错误
    WARNING = 2     # 疑似错误


class Detector:
    def __init__(self, cfg, wrkdir, bam_fname, asm_fname):
        self.cfg = cfg
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
        
        #merged = [("contig_5", 34400, 35000)]
        misassemblies = self.verify_candidates(merged)

        self.dump_misassemblies(os.path.join(self.wrkdir, "misassembly.bed"), misassemblies)


    @staticmethod
    def is_read_trivial(read, mapq):
        #return read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < mapq
        return read.is_unmapped or read.is_secondary or read.mapping_quality < mapq

    @staticmethod
    def detect_region_type(ctg, start, end, bam, ref):
        #check_win_size = self.cfg["check_win_size"]
        check_win_size = 5000
        #(lcov, rcov) = Detector.check_surrounding_region(ctg, start, end)
        #th_cov = max(3, min(lcov, rcov) / 3)
        th_cov = 3
        win_iter = utils.WinIterator(end - start, check_win_size, check_win_size)
        for s, e in win_iter:
            s += start
            e += start
            t = Detector.detect_win_type(ctg, s, e, bam, ref, th_cov)
            if t != RegionType.NORMAL:
                return t
        return RegionType.NORMAL
    
    @staticmethod
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
        

    @staticmethod
    def calc_distance(start, end, read, ref):
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

    @staticmethod
    def detect_win_type(ctg, start, end, bam, ref, th_cov):
        
        min_mapq = 10#self.cfg["min_mapq"]
        min_clip_len = 500#self.cfg["min_clip_len"]
        min_distance = 0.3#self.cfg["min_distance"]

        span = []
        clip_num = 0
        lowqual_num = 0
        cov = 0
        ctg_len = bam.get_reference_length(ctg)
        for read in bam.fetch(ctg, start, end):
            if (Detector.is_read_trivial(read, min_mapq)):
                continue
            cov += 1
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
            if Detector.calc_distance(start, end, read, ref) < 0.3:
                sstart, send = max(0, start - 500),  min(ctg_len, end + 500)
                if (read.reference_start <= sstart or sstart == 0 and read.reference_start < 500)  and \
                   (read.reference_end >= send or send == ctg_len and send + 500 > ctg_len):
                    span.append(read)
            else:
                lowqual_num += 1
        
        utils.logger.debug(f"detect_win_type: {ctg}:{start}-{end} {len(span)} >= {th_cov} {lowqual_num} {clip_num}")
        if len(span) < th_cov:
            return RegionType.ERROR
        elif clip_num >= cov // 3 or lowqual_num > cov // 3:
            utils.logger.debug(f"detect_win_type(WARNING): {ctg}:{start}-{end} {clip_num} >= {cov} {lowqual_num} > {cov}")
            return RegionType.WARNING
        else:
            return RegionType.NORMAL
    
    @staticmethod
    def verify_candidates_block(bam_fname, asm_fname, candidates, threads, i):
        verified = []
        utils.logger.info(f"start verify {i}/{threads}")

        bam = pysam.AlignmentFile(bam_fname, "rb")
        ref = pysam.FastaFile(asm_fname)
        bsize = (len(candidates) + threads - 1) // threads
        start = bsize * i
        end = min(len(candidates), bsize * (i+1))
        for ctg_name, s, e in candidates[start:end]:
            utils.logger.info(f"start verify: {ctg_name}:{s}-{e}")
            t = Detector.detect_region_type(ctg_name, s, e, bam, ref)
            if t != RegionType.NORMAL:
                verified.append((ctg_name, s, e, t))
        
        utils.logger.info("misassembly size = {}".format(len(verified)))
        return verified

    def verify_candidates(self, candidates):
        
        # threads_data = defaultdict(list)
        # for c in candidates:
        #     threads_data[c[0]].append([c. self.summary)

        threads = self.cfg["threads"]
        verified = []

        results = []
        pool = mp.Pool(processes=threads)
        for i in range(threads):
            results.append(pool.apply_async(Detector.verify_candidates_block, args=(self.bam_fname, self.asm_fname, candidates, threads, i)))

        pool.close() 
        pool.join()

        for r in results:
            print(r.get())
            verified.extend(r.get())
        
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
            for ctg_name, start, end, type in misassembies:
                f.write("%s\t%d\t%d\t%d\n" % (ctg_name, start, end, type.value))
