import os

from functools import partial
import multiprocessing as mp

import utils

from collections import Counter

class PileupInfo(object):
    
    def __init__(self, wrkdir):
        self.wrkdir = wrkdir
        self.pileups = {}

    def _get_pileup_fname(self, ctg):
        return os.path.join(self.wrkdir, "%s_pileup.txt" % ctg)
    
    def build(self, ctgs, bam_fname, asm_fname, threads, min_mapq):

        results = []
        pool = mp.Pool(processes=threads)
        for ctg_name, ctg_len in ctgs:
            results.append(pool.apply_async(partial(self._build_one_contig), args=(ctg_name, bam_fname, asm_fname, min_mapq)))
        pool.close() 
        pool.join()

        for r in results:
            ctg_name, pileup = r.get()
            self.pileups[ctg_name] = pileup


    def _build_one_contig(self, ctg_name, bam_fname, asm_fname, mapq):
        pileup_fname = self._get_pileup_fname(ctg_name)
        #utils.run_samtools_mpileup(bam_fname, pileup_fname, ctg_name, asm_fname, mapq)

        return ctg_name, self.load_pileup(pileup_fname)


    def load_pileup(self, fname):
        infos = []
        for line in open(fname):
            its = line.split()
            assert len(its) == 6

            depth = int(its[3])
            detail = its[4]
            counter = Counter(detail)
            match = counter['.'] + counter[',']
            diff = counter['A'] + counter['C']+counter['G'] + counter['T'] + \
                   counter['a'] + counter['c']+counter['g'] + counter['t']

            infos.append((depth, match, diff))
        return infos
    
    def get_mis_candidates(self, win_size, stride, max_diff_ratio):
        candidates = []
        for ctg_name, pileup in self.pileups.items():
            compressed = self.compress(pileup, win_size, stride)

            for c, (s, e) in zip(compressed, utils.WinIterator(len(pileup), win_size, stride)):
                 if c > max_diff_ratio:
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
                b += pileup[ii][2]
            
            compressed[i] = b / a if a > 0 else 1

        return compressed