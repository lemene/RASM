import os
import gzip
from collections import defaultdict
import itertools

import numpy as np

import utils
from feature import *

class DepthInfo(Feature):
    class Info(object):
        def __init__(self, chr_id, ctg_len):
            
            self.ctg_len = ctg_len
            self.chr_id = chr_id
            self.depths = []
            
    
    def __init__(self, wrkdir, win_size=100, mapq=10):
        self.wrkdir = wrkdir
        self.win_size = win_size
        self.mapq = mapq
        self.depths = {}
        self.lengths = {}
        self.median = 0.0

    def build(self, ctgs, bam, threads):
                
        #utils.run_mosdepth(bam, self.wrkdir + "/depth", self.win_size, self.mapq, threads)
        depth_fname = os.path.join(self.wrkdir, "depth.regions.bed.gz")
        assert os.path.exists(depth_fname) and "The 'xxx.regions.bed.gz' should be generated"

        for ctg_name, ctg_len in ctgs:
            self.depths[ctg_name] = []
            self.lengths[ctg_name] = ctg_len
            
        for line in gzip.open(depth_fname, "rt"):
            its = line.split()
            if its[0] in self.depths:
                self.depths[its[0]].append(float(its[3]))
        
        self.median = np.median(sum(self.depths.values(), []))
        utils.logger.info("Global depth: %.02f" % self.median)

    def get_mis_candidates(self, win_size, stride, min_radio, max_radio):
        candidates = []
        for ctg_name, depth in self.depths.items():

            compressed = self.compress(depth, win_size, stride)

            count = 0
            for d, (s, e) in zip(compressed, utils.WinIterator(self.lengths[ctg_name], win_size, stride)):
                count += 1
                #print(d, min_radio * self.median, max_radio * self.median)
                if d < min_radio * self.median or d > max_radio * self.median:
                    candidates.append((ctg_name, s, e))   

        for ctg_name, start, end in candidates:
            utils.logger.info("cand depth: %s:%d-%d" % (ctg_name, start, end))
        return candidates

    def compress(self, depth, win_size, stride):
        nwin_size = win_size // self.win_size
        nstride = stride // self.win_size

        win_iter = utils.WinIterator(len(depth), nwin_size, nstride)
        compressed = [0] * win_iter.size()

        for i, (s, e) in enumerate(win_iter):
            compressed[i] = sum(depth[s:e]) / (e - s)

        return compressed