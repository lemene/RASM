
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
        clips = []
        quals = []
        for ctg, (table, reads) in self.infos.items():
            for accu, lclip, rclip in reads.values():
                if lclip > 5: clips.append(lclip)
                if rclip > 5: clips.append(rclip)
                quals.append(accu)
        s = scipy.stats.median_abs_deviation(quals)

        quals.sort()
        m = np.median(quals)
        print("s,m,", s, m, quals[len(quals)//2])

        bin_clips = np.zeros(10000)
        for c in clips:
            if c < len(bin_clips):
                bin_clips[c] += 1

        for i in range(1000-100):
            print(i, sum(bin_clips[i:i+100])/100)


    def stat(self, threads):
        #self.stat_read()
        self.stat_table(threads)
        #self.stat_depth(threads)

    @staticmethod
    def stat_region_depth(table):

        depths = [0]*1000
        others = defaultdict(int)
        for t in table.T:
            s = int(sum(t[0:5]))
            if s < len (depths):
                depths[s] += 1
            else:
                others[s] += 1
        return depths, others
            
    def stat_depth(self, threads):
        utils.logger.info("Start stat depth")
        N = 10000000

        results = []
        pool = mp.Pool(processes=threads)
        for ctg, (table, reads) in self.infos.items():
            start = 0
            while start < table.shape[1]:
                end = min(start + N, table.shape[1])
                
                results.append(pool.apply_async(Summary.stat_region_depth, args=(table[:,start:end],)))
                start = end
        
        depths = [0]*1000
        others = defaultdict(int)
        for r in results:
            deps, oths = r.get()

            for i, d in enumerate(deps):
                depths[i] += d
            
            for k, v in oths.items():
                others[k] += v

        others = sorted(others.items(), key=lambda x: x[0])
        print("others", len(others))

        # calucate median
        count = sum(depths) + sum([i[1] for i in others])
        accu_m = 0
        median = 0
        for d, c in enumerate(depths):
            accu_m += c
            if accu_m*2 >= count:
                median = d
                break
        else:
            for d, c in others.item():
                accu_m += c
                if accu_m*2 >= count:
                    median = d
                    break

        utils.logger.info("median of depths: %d" % median)

        assert median < len(depths)
        K = 5
        smooths = [0] * (min(len(depths), median) - K + 1)
        smooths[0] = sum(depths[0:K])
        for i in range(1,len(smooths)):
            smooths[i] = smooths[i-1] - depths[i-1] + depths[i-1+K]

        mx = np.argmin(smooths)
        return median, mx

    #def get_mis_cands_from_depth(self, threads):



    def stat_table(self, threads):
        bsize = 10000000
        results = []
        pool = mp.Pool(processes=threads)
        for ctg, (table, reads) in self.infos.items():
 
            s = 0
            tablen = table.shape[1]
            while s < tablen:
                e = min(s + bsize, tablen)
                results.append(pool.apply_async(Summary.stat_region, args=(table[:,s:e],)))
                s = e

        pool.close() 
        pool.join()
        ref_total, ref_ins, ref_ok, rd_total, rd_ins, rd_ok = 0, 0, 0, 0, 0, 0
        for r in results:
            rr = r.get()
            ref_total += rr[0]
            ref_ins += rr[1]
            ref_ok += rr[2]
            rd_total += rr[3]
            rd_ins += rr[4]
            rd_ok += rr[5]
            print("sub",rr[2] / (rr[0] + rr[1]), rr[5] / (rr[3] + rr[4]))

        print("totoal", ref_ok / (ref_total + ref_ins), rd_ok / (rd_total + rd_ins))

    @staticmethod
    def stat_region(table):
        
        # accu = match / (match + mismatch + insertion + deletion)
        ref_ok, ref_err, ref_total, ref_ins = 0, 0, 0, 0
        rd_ok, rd_err, rd_total, rd_ins = 0, 0, 0, 0
        for i, t in enumerate(table.T):
            s = sum(t[0:8])
            mx = np.argmax(t[0:8])
            
            if mx == 0:
                ref_ok += 1
                rd_err += sum(t[1:5])
                rd_ok += t[0]
                #rd_total += sum(t[0:5])
            elif mx >=1 and mx <= 4:
                ref_err += 1
                rd_ok += t[mx]
                rd_err += sum(t[0:5]) - t[mx]
                #rd_total += sum(t[0:5])
            elif mx == 5:               # deletion
                #rd_err += sum(t[0:5])
                #assert 0
                rd_ins += sum(t[0:5])
                pass
            elif mx == 6:               # insertion
                #rd_err += sum(t[0:5])
                #assert 0
                ref_ins += t[7] / t[6]
                #rd_ok += t[7]
                pass

            rd_total += sum(t[0:5]) + t[7] 
            rd_ins += t[6]

        ref_total += len(table.T)
        
        utils.logger.info("done: ref_total=%d, ref_ins=%d, ref_ok=%d, rd_total=%d, rd_ins=%d, rd_ok=%d" % \
                           (ref_total, ref_ins, ref_ok, rd_total, rd_ins, rd_ok))
        return ref_total, ref_ins, ref_ok, rd_total, rd_ins, rd_ok
        

    def split_contig_tasks(self, ctgs, bsize=10000000):

        segs = []
        for ctg_name, ctg_len in ctgs:
            s = 0
            while s < ctg_len:
                e = min(s+bsize, ctg_len)
                segs.append((ctg_name, ctg_len, s, e))
                s = e
        return segs


    def scan(self, threads, min_contig):
        ctgs = self.get_contigs(min_contig)
        
        segs = self.split_contig_tasks(ctgs, 1000000)
        utils.logger.info("summary: split ctgs to %d parts", len(segs))
        # a = self.scan_segment("ctg10", 912736, 456530, 456540)
        # self.stat_region(a[4])
        # assert 0

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
                self.infos[ctg_name] = [np.zeros([9, ctg_len]), {}]
            self.infos[ctg_name][0][:, start:end] = table
            self.infos[ctg_name][1].update(reads)
        
        return ctg_name, table, reads

    def scan_segment(self, ctg_name, ctg_len, start, end):
        '''扫描bam文件获得基础信息'''
        utils.logger.info("start scanning %s:%d-%d" % (ctg_name, start, end))
        bam = pysam.AlignmentFile(self.bam_fname, "rb")
        ref = pysam.FastaFile(self.asm_fname)
        is_trivial = lambda r: r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 1

        reads = {}
        table = np.zeros([9, end-start])        # 0 match，1-4 mismatch, 5 deletion 6 insertion 7 insertionsize
        rseq = ref.fetch(ctg_name).upper()
        for read in bam.fetch(ctg_name, start, end):
            if is_trivial(read): continue
            rpos = read.reference_start
            if len(rseq) != ctg_len:
                utils.logger.info("%s:-%d-%d %d %d" % (ctg_name, start, end, len(rseq), ctg_len))
                assert 0
            qpos = 0
            qseq = read.query_sequence.upper()
            Bases = {'A':0, 'C':1,'G':2,'T':3}

            lclip = 0 if read.cigartuples[0][0] != 4 and read.cigartuples[0][0] != 5 else read.cigartuples[0][1]
            rclip = 0 if read.cigartuples[-1][0] != 4 and read.cigartuples[-1][0] != 5 else read.cigartuples[-1][1]
            match, mismatch, insertion, deletion = 0, 0, 0, 0
    
            for op, op_len in read.cigartuples:
                if op == 0:         # M
                    for i in range(op_len):
                        assert(qseq != None)
                        assert(rseq != None)
                        if rpos+i >= len(rseq) or qpos + i >= len(qseq):
                            utils.logger.info("%d: %d/%d %d/%d " % (i, rpos, len(rseq), qpos, len(qseq)))
                            assert 0
                        if rseq[rpos+i] == qseq[qpos+i]:
                            if rpos + i >= start and rpos + i < end:
                                table[0, rpos+i-start] += 1
                            match += 1
                        else:
                            if rpos + i >= start and rpos + i < end:
                                table[1+Bases[qseq[qpos+i]], rpos+i - start] += 1
                            mismatch += 1
                    
                    qpos += op_len
                    rpos += op_len
                    
                elif op == 1:   # I insertion
                    if rpos >= start and rpos < end:
                        table[6, rpos-start] += 1
                        table[7, rpos-start] += op_len
                    qpos += op_len
                    insertion += op_len
                elif op == 2:   # D deletion
                    if rpos >= start:
                        table[5, rpos-start:rpos-start+min(end-start, op_len)] += 1
                    rpos += op_len
                    deletion += op_len
                elif op == 7:   # =
                    if rpos >= start:
                        table[0, rpos-start:rpos-start+min(end-start,op_len)] += 1
                    qpos += op_len
                    rpos += op_len
                    match += op_len
                elif op == 8:   # X
                    for i in range(op_len):
                        if rpos + i >= start and rpos + i < end:
                            table[1+Bases[qseq[qpos+i]],rpos+i-start] += 1
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