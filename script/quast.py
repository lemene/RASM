
import re
import sys
import os
import copy
from collections import defaultdict
# quast_file = sys.argv[1]
# bed = sys.argv[2]

class Align_Info:
    def __init__(self, contig, start, end, ref, ref_start, ref_end, IDY, seg_id) -> None:
        self.contig = contig
        self.start = start
        self.end = end
        ##
        self.ref = ref
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.IDY = IDY
        ##
        self.seg_id = seg_id
        self.is_left = False
        self.is_right = False
        #
        self.align_orient = "+"
        if self.start > self.end:
            self.align_orient = "-"
            self.start, self.end = self.end, self.start
        #
    # def

class Mis_rec:
    def __init__(self, ctg, start, end, seg_id_ls, unalign_len, mis_info_ls) -> None:
        self.ctg = ctg
        self.start = start
        self.end = end
        self.seg_id_ls = seg_id_ls
        self.unalign_len = unalign_len
        self.mis_info_ls = mis_info_ls
        ##
        self.ref = None
        self.ref_start = None
        self.ref_end = None

    def Add_ref(self, ref1, start1, end1, ref2, start2, end2, is_bound):
        if is_bound:
            self.ref = ref1
            self.ref_start = min(start1, end1)
            self.ref_end = max(start1, end1)
        elif ref1 == ref2:
            self.ref = ref1
            start1, end1 = min(start1, end1), max(start1, end1)
            start2, end2 = min(start2, end2), max(start2, end2)
            if start1 > start2:
                start1, end1, start2, end2 = start2, end2, start1, end1
            self.ref_start = min(end1, start2)
            self.ref_end = max(end1, start2)
        else:
            self.ref = []
    @staticmethod
    def write_rec(mis_rec_ls, bed):
        with open(bed, "w") as fo:
            fo.write("#ctg\tstart\tend\tseg_ids\tunalign_len\tmis_info_ls\n")
            for mis_rec in mis_rec_ls:
                # mis_rec = Mis_rec()
                fo.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(mis_rec.ctg, str(mis_rec.start), str(mis_rec.end), ",".join('%s' %id for id in mis_rec.seg_id_ls), str(mis_rec.unalign_len), ",".join(mis_rec.mis_info_ls)))
    @staticmethod
    def read_rec(bed):
        dic = defaultdict(list)
        with open(bed, "r") as f:
            for line in f:
                if line.startswith("#"): continue
                ctg, start, end, seg_ids, unalign_len, mis_info_ls = line.strip().split("\t")
                # print(start, end)
                start, end = int(start), int(end)
                # seg_id_ls = seg_ids.split(",")
                # unalign_len = int(unalign_len)
                # mis_info_ls = mis_info_ls.split(",")
                # mis_rec = Mis_rec(ctg, start, end, seg_id_ls, unalign_len, mis_info_ls)
                # dic[ctg].append(mis_rec)
                dic[ctg].append([start, end])
        return dic

''' contig types
3 unalign types
1、correct_unaligned
2、mis_unaligned
3、unaligned

4、correct      无misassembly
5、misassembled 有
'''
'''
1、indels
2、local misassembly
3、relocation、translocation、inversion
'''


'''
fake: not a misassembly (possible transposable element)
indel: indel (<= 5bp)
indel: indel (> 5bp)
indel: stretch of mismatches

local misassembly
relocation, inconsistency = -10081
translocation
inversion
'''
def convert4(tsv_file, bed, bias=5):

    with open(tsv_file, "r") as fin, open(bed, "w") as fo:
        print("Start")
        line = fin.readline()
        # print(line)
        line = fin.readline()
        ctg_align_info_ls = []
        ctg_mis_info_ls = []
        mis_rec_ls = []
        seg_id = 0
        # bias = 10
        unalign_num = 0
        correct_num = 0
        while line:
            # print(line)
            if line.startswith("CONTIG"):
                ## 处理last contig info
                fields = line.strip().split("\t")
                contig_type = fields[3]
                ctg, ctg_len = fields[1:3]
                ctg_len = int(ctg_len)
                if len(ctg_mis_info_ls) == 0 and (contig_type == "mis_unaligned" or contig_type == "unaligned" or contig_type == "correct_unaligned"):
                    print("---------{} is {}, no mis ---------".format(ctg, contig_type))
                    unalign_num += 1
                elif len(ctg_mis_info_ls) == 0 and contig_type == "correct":
                    print("---------{} is {}---------".format(ctg, contig_type))
                    correct_num += 1
                else:   # misassembled contig、correct contig(has local mis)
                    if contig_type == "mis_unaligned" or contig_type == "unaligned" or contig_type == "correct_unaligned":
                        print("{}, {}:{}".format(ctg, contig_type, ctg_mis_info_ls))
                    print("{} start, type:{}".format(ctg, contig_type))
                    ctg_align_info_ls[0].is_left = True
                    ctg_align_info_ls[-1].is_right = True
                    #
                    pre_align_info = None
                    mis_idx = 0
                    for align_info in ctg_align_info_ls:
                        # align_info = Align_Info()
                        if align_info.is_left:  # Type1: left
                            mis_rec = Mis_rec(ctg, 1, align_info.start, [align_info.seg_id], align_info.start - 1, ['left_unalign'])
                            # mis_rec.
                            # print("left unalign:1-{}, length:{}".format(str(align_info.start), str(align_info.start - 1)))
                            # print(mis_rec.__dict__)
                            mis_rec_ls.append(mis_rec)
                        elif ctg_mis_info_ls[mis_idx][0] == "unknown":  # unknown type
                            mis_idx += 1
                            continue
                        else:
                            # Type2
                            mis_info = ctg_mis_info_ls[mis_idx]
                            # pre_align_info = Align_Info()
                            # align_info = Align_Info()
                            mis_start = max(pre_align_info.end - bias, pre_align_info.start)
                            mis_end = min(align_info.start + bias, align_info.end)
                            seg_id_ls = [pre_align_info.seg_id, align_info.seg_id]
                            unalign_len = align_info.start - pre_align_info.end
                            mis_rec = Mis_rec(ctg, mis_start, mis_end, seg_id_ls, unalign_len, mis_info)
                            if mis_rec.start > mis_rec.end:
                                print("pre:{}-{}, now:{}-{}".format(pre_align_info.start, pre_align_info.end, align_info.start, align_info.end))
                                print("Error")
                                print(pre_align_info.__dict__, align_info.__dict__)
                                raise ValueError
                            mis_idx += 1
                            mis_rec_ls.append(mis_rec)
                            # Type3
                            if align_info.is_right: # 右边界的末端-end
                                mis_rec = Mis_rec(ctg, align_info.end, ctg_len, [align_info.seg_id], ctg_len - align_info.end, ['right_unalign'])
                                # print("right unalign:{}-{}, length:{}".format(str(align_info.end), ctg_len, str(ctg_len - align_info.end)))
                                mis_rec_ls.append(mis_rec)
                        pre_align_info = align_info
                #
                print(ctg + ":", ctg_len, "finished\n")
                # break
                seg_id = 0
                ctg_align_info_ls = []
                ctg_mis_info_ls = []
                line = fin.readline()
            else:
                if line.strip() == "unknown":
                    line = fin.readline()
                    continue
                ## parse alignsegment info
                align_seg_fields = line.strip().split("\t")
                # print(align_seg_fields)
                # 181', '12634', '140', '12550', 'chrI', 'NC_001133.9', '96.82', '', 'True
                # S1      E1      S2      E2      Reference       Contig  IDY     Ambiguous       Best_group
                # print(align_seg_fields)
                ref_start, ref_end, ctg_start, ctg_end, ref, contig, IDY = align_seg_fields[:7]
                # if int(ctg_start) > int(ctg_end):
                #     print(ref_start, ref_end, ctg_start, ctg_end, ref, contig, "is reverse")
                # print(ref_start, ref_end, ctg_start, ctg_end, ref, contig, IDY, seg_id)
                # contig, start, end, ref, ref_start, ref_end, IDY, seg_id
                align_info = Align_Info(contig, int(ctg_start), int(ctg_end), ref, int(ref_start), int(ref_end), float(IDY), seg_id)
                # print(align_info.__dict__)
                ctg_align_info_ls.append(align_info)
                line = fin.readline()
                if line.startswith("CONTIG"):   # End
                    continue
                else:   # parse mis info
                    mis_info = line.strip().split(",")
                    # mis_info = line.strip()
                    # print(mis_info)
                    ctg_mis_info_ls.append(mis_info)
                    line = fin.readline()
                seg_id += 1

        ## filter。主要是
        new_mis_rec_ls = []
        for mis_rec in mis_rec_ls:

            if mis_rec.mis_info_ls[0] == 'left_unalign' or mis_rec.mis_info_ls[0] == 'right_unalign':
                # if mis_rec.unalign_len > 10000:
                #     new_mis_rec_ls.append(mis_rec)
                continue
            elif mis_rec.mis_info_ls[0].startswith("fake") or mis_rec.mis_info_ls[0].startswith("indel"):
                continue
            else:
                new_mis_rec_ls.append(mis_rec)
        print("-----------Sum:---------\nmis_num:", len(new_mis_rec_ls))
        print("unalign contig:{}".format(unalign_num))
        print("correct contig:{}".format(correct_num))
        Mis_rec.write_rec(new_mis_rec_ls, bed)

if __name__ == "__main__":
    tsv_file = sys.argv[1]
    bed = sys.argv[2]
    bias = int(sys.argv[3])
    convert4(tsv_file, bed, bias)

                                                                                                                                                                                     