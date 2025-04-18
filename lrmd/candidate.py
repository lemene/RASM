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

    def get_mis_candidates(self, ctgs, smry, cfg):
        results = []
        pool = mp.Pool(processes=cfg["threads"])
        for binfo in utils.split_contig_by_block(ctgs):
            utils.logger.info(f"{binfo}")
            ctg_name, ctg_len, start, end = binfo
            results.append(pool.apply_async(Candidate.collect_in_block, 
                                            args=(binfo, smry.infos[ctg_name][0][start:end,:], smry.stats, cfg)))
        pool.close()
        pool.join()

        cands = []
        for r in results:
            cands.extend(r.get())
        return cands

        
    @staticmethod
    def collect_in_block(binfo, table, stats, cfg):
        win_size = cfg["win_size"]
        stride = cfg["stride"]
        max_diff_ratio = cfg["max_diff_ratio"]
        min_clip_num = cfg["min_clip_num"]

        ctg_name, ctg_len, start, end = binfo
        win_iter = utils.WinIterator(end-start, win_size, stride)
        
        cands = []
        for i, (s, e) in enumerate(win_iter):
            
            # 该区域的最小覆盖低
            if np.min(np.sum(table[s:e, 0:5],axis=1)) < 3:#
                cands.append([ctg_name, s + start, e + start])
                continue
        
            # 区域的相似度
            bs = np.sum(np.sum(table[s:e, 0:5],axis=1) + table[s:e, 5] + table[s:e, 7])
            ms = np.sum(table[s:e, 0])
            if ms / bs < 1 - max_diff_ratio:
                cands.append([ctg_name, s + start, e + start])
                continue

            clips = np.sum(table[s:e, 8])
            if clips > min_clip_num:
                cands.append([ctg_name, s + start, e + start])
                continue

        return cands