from collections import defaultdict
import numpy as np
import multiprocessing as mp

import utils
from feature import *

class DepthInfo(Feature):
            
    def __init__(self):
        self.depths = {}
        self.median = 0.0
        self.trough = 0

    def build(self, ctgs, smry, threads):

        results = []
        pool = mp.Pool(processes=threads)
        for binfo in utils.split_contig_by_block(ctgs):
            (ctg_name, ctg_len, start, end) = binfo
            results.append(pool.apply_async(DepthInfo.collect_depth_in_block, 
                                            args=(binfo, smry.infos[ctg_name][0][:,start:end])))
        pool.close() 
        pool.join()

        for ctg_name, ctg_len in ctgs:
            self.depths[ctg_name] = np.zeros(ctg_len)

        for r in results:
            (ctg_name, ctg_len, start, end), depths = r.get()
            self.depths[ctg_name][start:end] = depths

        self.median, self.trough = self.calc_median()
        utils.logger.info("Global depth: %.02f" % self.median)

    @staticmethod
    def collect_depth_in_block(binfo, block):
        (ctg_name, ctg_len, start, end) = binfo
        infos = np.zeros(end-start)

        for i, t in enumerate(block.T):
            infos[i] = sum(t[0:5])

        return binfo, infos

    def calc_median(self):
        hist_1000 = np.zeros(1000)
        hist_others = defaultdict(int)
        for ds in self.depths.values():
            for d in ds:
                if d < 1000:
                    hist_1000[int(d)] += 1
                else:
                    hist_others[int(d)] += 1
        
        hist_others = sorted(hist_others.items(), key=lambda x: x[0])

        # calucate median
        count = sum(hist_1000) + sum([i[1] for i in hist_others])
        accu_m = 0
        median = 0
        for d, c in enumerate(hist_1000):
            accu_m += c
            if accu_m*2 >= count:
                median = d
                break
        else:
            for d, c in hist_others.item():
                accu_m += c
                if accu_m*2 >= count:
                    median = d
                    break

        utils.logger.info("median of depths: %d" % median)

        assert median < len(hist_1000)
        K = 5
        smooths = [0] * (min(len(hist_1000), median) - K + 1)
        smooths[0] = sum(hist_1000[0:K])
        for i in range(1,len(smooths)):
            smooths[i] = smooths[i-1] - hist_1000[i-1] + hist_1000[i-1+K]

        trough = np.argmin(smooths)
        utils.logger.info("first trough of depths: %d" % trough)
        return median, trough

    def get_mis_candidates(self, win_size, stride):
        candidates = []
        for ctg_name, depth in self.depths.items():

            compressed = self.compress(depth, win_size, stride)

            for d, (s, e) in zip(compressed, utils.WinIterator(len(depth), win_size, stride)):
                if d < self.trough:
                    candidates.append((ctg_name, s, e))   

        for ctg_name, start, end in candidates:
            utils.logger.debug("cand depth: %s:%d-%d" % (ctg_name, start, end))
        return candidates

    def compress(self, depth, win_size, stride):
        win_iter = utils.WinIterator(len(depth), win_size, stride)
        compressed = [0] * win_iter.size()

        for i, (s, e) in enumerate(win_iter):
            compressed[i] = min(depth[s:e])
        return compressed