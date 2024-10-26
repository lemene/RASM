import pysam
from functools import partial
import multiprocessing as mp

import utils
from feature import *

class ClipInfo(Feature):
    def __init__(self):
          self.clips = {}

    def build(self, ctgs, bam_fname, threads, min_mapq, min_clip_len):
        
        results = []
        pool = mp.Pool(processes=threads)
        for ctg_name, ctg_len in ctgs:
            results.append(pool.apply_async(partial(self._build_one_contig), args=(ctg_name, bam_fname, min_mapq, min_clip_len)))
        pool.close() 
        pool.join()

        for r in results:
            ctg_name, clip = r.get()
            self.clips[ctg_name] = clip


    def _build_one_contig(self, ctg_name, bam_fname, min_mapq, min_clip_len):
        bam = pysam.AlignmentFile(bam_fname, "rb")

        ctg_len = bam.get_reference_length(ctg_name)
        clips = [0] * ctg_len
        for read in bam.fetch(ctg_name):
            if read.is_unmapped or read.is_secondary or read.mapping_quality < min_mapq:
                continue

            left = read.cigartuples[0]
            if (left[0] == 4 or left[0] == 5)  and left[1] >= min_clip_len:
                    clips[read.reference_start] += 1

            right = read.cigartuples[-1]
            if (right[0] == 4 or right[0] == 5) and right[1] >= min_clip_len:
                    clips[read.reference_end-1] += 1

        return ctg_name, clips
         

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
        compressed = [0] * win_iter.size()

        for i, (s, e) in enumerate(win_iter):
            compressed[i] = sum(clips[s:e])

        return compressed