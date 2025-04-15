
import pysam
import numpy as np
import multiprocessing as mp
import functools as ft
import pickle
import scipy
from collections import defaultdict

import utils

class Summary(object):
    def __init__(self, bam_fname, asm_fname):
        self.bam_fname = bam_fname
        self.asm_fname = asm_fname
        self.infos = {}
        self.asm_accuracy = 0.0
        self.read_accuarcy = 0.0
        self.stats = {}

        self.depth_median = 0
        self.depth_trough = 0

    def load(self, fname):
        utils.logger.info("Load intermediate state: %s" % fname)
        self.infos = pickle.load(utils.open_file(fname, "rb"))

    def save(self, fname):
        with utils.open_file(fname, 'wb') as file:
            pickle.dump(self.infos, file)

    def get_contigs(self, min_contig):
        bam = pysam.AlignmentFile(self.bam_fname, "rb")
        ctgs = [(ref, rlen) for ref, rlen in zip(bam.references, bam.lengths) if rlen >= min_contig]
        return ctgs
    
    def stat_read(self):
        quals = []
        for ctg, (table, reads) in self.infos.items():
            for accu, lclip, rclip in reads.values():
                quals.append(accu)
        s = scipy.stats.median_abs_deviation(quals)

        quals.sort()
        m = np.median(quals)
        self.stats["read_qual_median"] = np.median(quals)
        self.stats["read_qual_max"] = quals[-1]
        self.stats["read_qual_min"] = quals[0]
        self.stats["read_qual_mad"] = s


    def stat(self, threads):
        utils.logger.info("stat_read")

        self.stat_read()
        utils.logger.info("stat_table")
        self.stat_table(threads)
        utils.logger.info("stat_depth")
        self.stat_depth(threads)

        for k in sorted(self.stats.values()):
            print(k, self.stats[k])
            utils.logger.info(f"{k} = {self.stats[k]}")

    def calc_depth_median(self, hist_1000, hist_gt1000):
        # calucate median
        count = sum(hist_1000) + sum([i[1] for i in hist_gt1000])

        accu_m = 0
        median = 0
        for d, c in enumerate(hist_1000):
            accu_m += c
            if accu_m*2 >= count:
                median = d
                break
        else:
            for d, c in hist_gt1000.item():
                accu_m += c
                if accu_m*2 >= count:
                    median = d
                    break

        return median

    def calc_depth_trough(self, hist):
        K = 5
        smooths = [0] * (len(hist) - K + 1)
        smooths[0] = sum(hist[0:K])
        for i in range(1,len(smooths)):
            smooths[i] = smooths[i-1] - hist[i-1] + hist[i-1+K]

        return np.argmin(smooths)
            
    def stat_depth(self, threads):
        utils.logger.info("Start stat depth")

        hist_1000 = [0]*1000
        hist_gt1000 = defaultdict(int)
        for ctg, (table, reads) in self.infos.items():
            for t in table:
                s = int(sum(t[0:5]))
                if s < len(hist_1000):
                    hist_1000[s] += 1
                else:
                    hist_gt1000[s] += 1

        hist_gt1000 = sorted(hist_gt1000.items(), key=lambda x: x[0])


        median = self.calc_depth_median(hist_1000, hist_gt1000)
        self.stats["depth_median"] = median

        assert median < len(hist_1000)
        trough = self.calc_depth_trough(hist_1000[0:median])
        self.stats["depth_trough"] = trough


    def stat_table(self, threads):

        ref_match, ref_total = 0, 0
        rd_match, rd_total = 0, 0
        for ctg, (table, reads) in self.infos.items():

            for _, t in enumerate(table):

                s = sum(t[0:8])
                mx = np.argmax(t[0:7])
                assert mx != 7
            
                if mx == 0:
                    ref_match += 1
                    rd_match += t[0]
                elif mx >=1 and mx <= 4:
                    rd_match += t[mx]
                elif mx == 5:               # deletion
                    pass
                elif mx == 6:               # insertion
                    rd_match += t[7]

                rd_total += sum(t[0:5]) + t[7] 

            ref_total += len(table)
        self.stats["ref_accu"] = ref_match / ref_total
        self.stats["read_accu"] = ref_match / rd_total


    def scan(self, threads, min_contig):
        ctgs = self.get_contigs(min_contig)
        
        segs = utils.split_contig_by_block(ctgs)

        results = []
        pool = mp.Pool(processes=threads)
        for ctg_name, ctg_len, start, end in segs:
            results.append(pool.apply_async(ft.partial(self.scan_segment), args=(ctg_name, ctg_len, start, end)))

        pool.close() 
        pool.join()

        utils.logger.info("Merge scanning results")
        for r in results:
            ctg_name, ctg_len, start, end, table, reads = r.get()
            if ctg_name not in self.infos:
                self.infos[ctg_name] = [np.zeros([ctg_len,9]), {}]
            self.infos[ctg_name][0][start:end,:] = table
            self.infos[ctg_name][1].update(reads)

        self.stat(threads)

    def scan_segment(self, ctg_name, ctg_len, start, end):
        '''扫描bam文件获得基础信息'''
        utils.logger.info("start scanning %s:%d-%d" % (ctg_name, start, end))
        bam = pysam.AlignmentFile(self.bam_fname, "rb")
        ref = pysam.FastaFile(self.asm_fname)
        is_trivial = lambda r: r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 0

        reads = {}
        table = np.zeros([end-start, 9])        # 0 match，1-4 mismatch, 5 deletion 6 insertion 7 insertionsize
        rseq = ref.fetch(ctg_name).upper()
        
        for read in bam.fetch(ctg_name, start, end):
            if is_trivial(read): continue
            rpos = read.reference_start
            qpos = 0
            qseq = read.query_sequence.upper()
            Bases = {'A':0, 'C':1,'G':2,'T':3}

            lclip = 0 if read.cigartuples[0][0] != 4 and read.cigartuples[0][0] != 5 else read.cigartuples[0][1]
            rclip = 0 if read.cigartuples[-1][0] != 4 and read.cigartuples[-1][0] != 5 else read.cigartuples[-1][1]
            match, mismatch, insertion, deletion = 0, 0, 0, 0
    
            for op, op_len in read.cigartuples:
                if op == 0:         # M
                    for i in range(op_len):

                        if rseq[rpos+i] == qseq[qpos+i]:
                            if rpos + i >= start and rpos + i < end:
                                table[rpos+i-start, 0] += 1
                            match += 1
                        else:
                            if rpos + i >= start and rpos + i < end:
                                table[rpos+i - start, 1+Bases[qseq[qpos+i]]] += 1
                            mismatch += 1
                    
                    qpos += op_len
                    rpos += op_len
                    
                elif op == 1:   # I insertion
                    if rpos >= start and rpos < end:
                        table[rpos-start,6] += 1
                        table[rpos-start,7] += op_len
                    qpos += op_len
                    insertion += op_len
                elif op == 2:   # D deletion
                    if rpos >= start:
                        table[rpos-start:rpos-start+min(end-start, op_len),5] += 1
                    rpos += op_len
                    deletion += op_len
                elif op == 7:   # =
                    if rpos >= start:
                        table[rpos-start:rpos-start+min(end-start,op_len),0] += 1
                    qpos += op_len
                    rpos += op_len
                    match += op_len
                elif op == 8:   # X
                    for i in range(op_len):
                        if rpos + i >= start and rpos + i < end:
                            table[rpos+i-start, 1+Bases[qseq[qpos+i]]] += 1
                    qpos += op_len
                    rpos += op_len
                    mismatch += op_len
                elif op == 4:   # S
                    qpos += op_len
                elif op == 5:   # H
                    pass
                else:
                    pass        # 不处理
            accu = match / (read.query_length)
            reads[read.query_name] = [accu, lclip, rclip]
        utils.logger.info("end scanning %s:%d-%d" % (ctg_name, start, end))
        return ctg_name, ctg_len, start, end, table, reads
    
    def get(self, ctg, start, end):
        return self.infos[ctg][0][start:end,:]