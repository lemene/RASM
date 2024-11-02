import pysam
from functools import partial
import multiprocessing as mp

import numpy as np

import utils
from feature import *

class ClipInfo(Feature):
    def __init__(self):
          self.clips = {}

    def build(self, ctgs, bam_fname:str, threads:int, min_clip_len:int):
        
        results = []
        pool = mp.Pool(processes=threads)
        for binfo in utils.split_contig_by_block(ctgs):
            results.append(pool.apply_async(ClipInfo.collect_clip_in_block, args=(binfo, bam_fname, min_clip_len)))
        pool.close() 
        pool.join()

        for ctg_name, ctg_len in ctgs:
            self.clips[ctg_name] = np.zeros(ctg_len)

        for r in results:
            (ctg_name, ctg_len, start, end), clip = r.get()
            self.clips[ctg_name][start:end] = clip


    @staticmethod
    def collect_clip_in_block(binfo, bam_fname, min_clip_len):
        ctg_name, ctg_len, start, end = binfo
        
        bam = pysam.AlignmentFile(bam_fname, "rb")

        clips =  np.zeros(end - start)
        for read in bam.fetch(ctg_name, start, end):
            if read.is_unmapped or read.is_secondary or read.mapping_quality < 1:
                continue

            left = read.cigartuples[0]
            if (left[0] == 4 or left[0] == 5)  and left[1] >= min_clip_len and \
                read.reference_start >= start and read.reference_start < end:
                
                clips[read.reference_start - start] += 1

            right = read.cigartuples[-1]
            if (right[0] == 4 or right[0] == 5) and right[1] >= min_clip_len and \
                read.reference_end-1 >= start and read.reference_end-1 < end:

                clips[read.reference_end-1 - start] += 1

        return binfo, clips
         

    def get_mis_candidates(self, win_size, stride, min_clip_num):
        candidates = []
        for ctg_name, clips in self.clips.items():
            compressed = self.compress(clips, win_size, stride)

            for c, (s, e) in zip(compressed, utils.WinIterator(len(clips), win_size, stride)):
                 if c > min_clip_num:
                      candidates.append((ctg_name, s, e))
                      
        for ctg_name, start, end in candidates:
            utils.logger.info("cand clip: %s:%d-%d" % (ctg_name, start, end))
        return candidates

    def compress(self, clips, win_size, stride):
        win_iter = utils.WinIterator(len(clips), win_size, stride)
        compressed = np.zeros(win_iter.size())

        for i, (s, e) in enumerate(win_iter):
            compressed[i] = sum(clips[s:e])

        return compressed