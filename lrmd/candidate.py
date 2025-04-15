from collections import defaultdict
import numpy as np
import multiprocessing as mp

import utils
from feature import *

class Candidate(object):
            
    def __init__(self):
        self.depths = {}
        self.median = 0.0
        self.trough = 0

    def get_mis_candidates(self, ctgs, smry, win_size, stride, threads):
        results = []
        pool = mp.Pool(processes=threads)
        for binfo in utils.split_contig_by_block(ctgs):
            ctg_name, ctg_len, start, end = binfo
            results.append(pool.apply_async(Candidate.collect_in_block, 
                                            args=(binfo, smry.infos[ctg_name][0][:,start:end], win_size, stride)))
        pool.close()
        pool.join()

        
    @staticmethod
    def collect_in_block(binfo, smry, win_size, stride):
        ctg_name, ctg_len, start, end = binfo
        win_iter = utils.WinIterator(end-start, win_size, stride)
        
        cands = []
        for i, (s, e) in enumerate(win_iter):
            
        return cands