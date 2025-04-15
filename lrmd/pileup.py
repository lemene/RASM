import os

from functools import partial
import multiprocessing as mp

import utils
import numpy as np

from collections import Counter

class PileupInfo(object):
    
    def __init__(self):
        self.pileups = {} 

    def _get_pileup_fname(self, ctg):
        return os.path.join(self.wrkdir, "%s_pileup.txt" % ctg)
    
    def build(self, ctgs, smry, threads):

        results = []
        pool = mp.Pool(processes=threads)
        for binfo in utils.split_contig_by_block(ctgs):
            (ctg_name, ctg_len, start, end) = binfo
            
            results.append(pool.apply_async(PileupInfo.collect_pileup_in_block, 
                                            args=(binfo, smry.infos[ctg_name][0][:,start:end])))
        pool.close() 
        pool.join()

        
        for ctg_name, ctg_len in ctgs:
            self.pileups[ctg_name] = np.zeros((ctg_len,3))

        for r in results:
            (ctg_name, ctg_len, start, end), pileup = r.get()
            print("pileup: ", start, end, pileup)
            self.pileups[ctg_name][start:end,:] = pileup

    @staticmethod
    def collect_pileup_in_block(binfo, block):
        (ctg_name, ctg_len, start, end) = binfo
        infos = np.zeros((end-start, 3))

        for i, t in enumerate(block.T):
            infos[i,0] = sum(t[0:5]) + t[5] + t[7]
            infos[i,1] = t[0]

        return binfo, infos
    
    def get_mis_candidates(self, win_size, stride, max_diff_ratio):
        candidates = []
        for ctg_name, pileup in self.pileups.items():
            compressed = self.compress(pileup, win_size, stride)

            for c, (s, e) in zip(compressed, utils.WinIterator(len(pileup), win_size, stride)):
                 
                print("check", s, e, c, max_diff_ratio)
                if 1-c > max_diff_ratio:
                    candidates.append((ctg_name, s, e)) 

        for ctg_name, start, end in candidates:
            utils.logger.info("cand pileup: %s:%d-%d" % (ctg_name, start, end))
        return candidates
            
    def compress(self, pileup, win_size, stride):
        win_iter = utils.WinIterator(len(pileup), win_size, stride)
        compressed = [0] * win_iter.size()

        for i, (s, e) in enumerate(win_iter):
            a, b = 0, 0
            for ii in range(s, e):
                a += pileup[ii][0]
                b += pileup[ii][1]
            
            compressed[i] = b / a if a > 0 else 1

        return compressed