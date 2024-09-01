
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
def merge_reg(ls):
    new_ls = []
    l, r = -1, -1
    for reg in ls:
        if l < 0:
            l, r = reg
        else:
            if reg[0] <= r:
                r = reg[1]
            else:
                new_ls.append([l, r])
                l, r = reg
    if l > 0:
        new_ls.append([l, r])
    return new_ls

def convert(quast_file, bed):
    print("Start convert")
    misassembly_pattern = re.compile("Extensive misassembly \((.*)\).*between ([0-9]*) ([0-9]*) and ([0-9]*) ([0-9]*)")
    with open(quast_file, "r") as fin, open(bed, "w") as fou:
        lines = fin.readlines()
        for line in lines:
            if not line: continue
            match = misassembly_pattern.match(line)
            # print(match)
            if(match is None):
                contig_name = line.split("\n")
                current_contig = contig_name[0]
            else:
                misassembly_type = match.group(1)
                #Assign coordinate values and swap them if they are in an ascending order
                mis_contig_1a = match.group(2)
                mis_contig_1b = match.group(3)
                if (int(mis_contig_1a) > int(mis_contig_1b)):
                    mis_contig_1a, mis_contig_1b = mis_contig_1b, mis_contig_1a
                    first_piece_inverted = True
                mis_contig_2a = match.group(4)
                mis_contig_2b = match.group(5)
                if (int(mis_contig_2a) > int(mis_contig_2b)):
                    mis_contig_2a, mis_contig_2b = mis_contig_2b, mis_contig_2a
                    second_piece_inverted = True
                # print(current_contig, mis_contig_1a, mis_contig_1b, mis_contig_2a, mis_contig_2b)
                # print("misassembly_type: ", misassembly_type)
                mis_contig_1a, mis_contig_1b, mis_contig_2a, mis_contig_2b = int(mis_contig_1a), int(mis_contig_1b), int(mis_contig_2a), int(mis_contig_2b)
                mis_l = max(mis_contig_1b - 1000, 1)
                mis_r = min(mis_contig_2b, mis_contig_2a + 1000)
                print(mis_l, mis_r)
                fou.write("{}\t{}\t{}\n".format(current_contig, str(mis_l), str(mis_r)))


def convert2(html_file, bed):
    with open(html_file, "r") as fin, open(bed, "w") as fo:
        chr_len = {}
        chr_ls = []
        bias = 50
        chr_mis_ls = []
        first = True
        now_left, now_right = 1, 1
        pre_dic = {}
        for line in fin:
            if line is None: continue
            if line.startswith("chromosomes_len"):
                input_string = line.strip("\n").strip(";")
                # chromosomes_len["NC_001133.9"] = 230218
                # references_by_id["0"] = "NC_001133.9"
                regex = r'chromosomes_len\["([^"]+)"\] = (\d+)'     
                match = re.match(regex, input_string)
                # print(match.group(1), match.group(2))
                chr = match.group(1).strip("\"")
                length = int(match.group(2))
                chr_len[chr] = length
                chr_ls.append(chr)
            if not line.startswith("{name"): continue
            if first:
                idx = 0
                now_chr = chr_ls[idx]
                now_chr_len = chr_len[now_chr]
                now_left = now_right
                now_right = now_left + now_chr_len
                first = False
                pre_chr = None
                print("Start", now_chr)
            
            # if line.startswith("{name"): print()
            fields = line.strip().strip(",").lstrip("{").rstrip("}").split(",")
            dic = {}
            for i in fields:
                key, val = i.split(":")
                if key == "misassemblies" or key == "mis_ends":
                    val = val.strip("\"")
                    if len(val) > 0:
                        val_ls = val.strip("\"").strip("\"").split(";")
                        # print(dic["name"], dic["start"], dic["end"], val, val_ls, len(val))
                    else:
                        val_ls = []
                        # print(dic["name"], dic["start"], dic["end"], val, val_ls, len(val))
                    dic[key] = val_ls
                elif key == "name" or key == "is_best":
                    val = val.strip("\"")
                    dic[key] = val
                else:
                    dic[key] = val
            
            if int(dic["corr_start"]) >= now_right or int(dic["corr_end"]) >= now_right:
                '''将上个染色体的存起来'''
                if len(chr_mis_ls) > 0:
                    print("Write", now_chr)
                    for mis in chr_mis_ls:
                        # print(now_chr, str(mis[0]), str(mis[1]))
                        fo.write("{}\t{}\t{}\n".format(now_chr, str(mis[0]), str(mis[1])))
                chr_mis_ls = []
                ##
                idx += 1
                now_chr = chr_ls[idx]
                now_chr_len = chr_len[now_chr]
                now_left = now_right
                now_right = now_left + now_chr_len
                print("Start", now_chr)
                
            if len(dic["misassemblies"]) > 0:
                mis_ls = dic["misassemblies"]
                mis_ends_ls = dic["mis_ends"]
                if pre_chr == now_chr:
                    pre_mis_ends_ls = pre_dic["mis_ends"]
                    '''if 'L' in mis_ends_ls:
                        mis_l = max(int(dic["start"]) - bias, 0)
                        mis_r = min(int(dic["start"]) + bias, now_chr_len)
                        chr_mis_ls.append([mis_l, mis_r])
                    if 'R' in mis_ends_ls:
                        mis_l = max(int(dic["end"]) - bias, 0)
                        mis_r = min(int(dic["end"]) + bias, now_chr_len)
                        chr_mis_ls.append([mis_l, mis_r])'''
                    if 'L' in mis_ends_ls and 'R' in pre_mis_ends_ls:
                        # mis_l = max(int(pre_dic["end"]) - bias, 0)
                        # mis_r = min(int(dic["start"]) + bias, now_chr_len)
                        l, r = min(int(pre_dic["end"]), int(dic["start"])), max((int(pre_dic["end"]), int(dic["start"])))
                        mis_l = max(l - bias, 0)
                        mis_r = min(r + bias, now_chr_len)
                        chr_mis_ls.append([mis_l, mis_r])
                    elif 'L' in mis_ends_ls:
                        print(now_chr, dic["start"], dic["end"], "only L now")
                    elif 'R' in pre_mis_ends_ls:
                        print(now_chr, dic["start"], dic["end"], "only R pre")
                else:
                    # print(now_chr, pre_chr)
                    if 'L' in mis_ends_ls:
                        mis_l = max(int(dic["start"]) - bias, 0)
                        mis_r = min(int(dic["start"]) + bias, now_chr_len)
                        chr_mis_ls.append([mis_l, mis_r])
                    if len(pre_dic) > 0:
                        pre_mis_ends_ls = pre_dic["mis_ends"]
                        if 'R' in pre_mis_ends_ls:
                            mis_l = int(pre_dic["end"]) - bias
                            mis_r = chr_len[pre_chr]
                            # print(("{}\t{}\t{}\n".format(pre_chr, mis_l, mis_r)))
                            fo.write("{}\t{}\t{}\n".format(pre_chr, mis_l, mis_r))
            
            pre_dic = dic
            pre_chr = now_chr


            
            # print(dic)
            # if int(dic["start"]) > int(dic["end"]):print(dic["name"], dic["start"], dic["end"])
            # pre_dic = dic
            # print(pre_dic)
            # print(pre_dic["name"])
            '''if not pre_dic:
                if int(dic["start"]) > 1:
                    mis_ctg = dic["name"].strip("\"")
                    # mis_l = dic["start"]
                    mis_l = int(dic["start"]) - bias
                    mis_r = int(dic["start"]) + bias
                    fo.write("{}\t{}\t{}\n".format(mis_ctg, str(mis_l), str(mis_r)))
                pre_dic = dic
                continue
            if dic["name"] != pre_dic["name"]:
                if int(dic["start"]) > 1:
                    mis_ctg = dic["name"].strip("\"")
                    # mis_l = dic["start"]
                    mis_l = int(dic["start"]) - bias
                    mis_r = int(dic["start"]) + bias
                    fo.write("{}\t{}\t{}\n".format(mis_ctg, str(mis_l), str(mis_r)))
                pre_dic = dic
                continue
            mis_ctg = dic["name"].strip("\"")
            # mis_l = dic["start"]
            mis_l = max(int(pre_dic["end"]) - bias, int(pre_dic["start"]))
            mis_r = min(int(dic["start"]) + bias, int(dic["end"]))
            fo.write("{}\t{}\t{}\n".format(mis_ctg, str(mis_l), str(mis_r)))
            pre_dic = dic'''
        if len(chr_mis_ls) > 0:
            print("Write", now_chr)
            for mis in chr_mis_ls:
                # print(now_chr, str(mis[0]), str(mis[1]))
                fo.write("{}\t{}\t{}\n".format(now_chr, str(mis[0]), str(mis[1])))
        print(chr_len)

def convert3(html_file, bed):   # parse contig_size_viewer.html
    print("Start")
    with open(html_file, "r") as fin, open(bed, "w") as fo:
        first = True
        line = fin.readline()
        ref_dic = {}
        while line:
            # print(line)
            if line.startswith("contig_structures") and first:
                first = False
                line = fin.readline() 
                line = fin.readline()
            elif line.startswith("contig_structures"):  ## 记录按照contig来存的，记载了某条contig的比对结构信息
                # parse for asm contig name
                input_string = line.strip()
                regex = r'contig_structures\["([^"]+)"\]\["([^"]+)"\] = \['
                match = re.match(regex, input_string)
                ref_id = match.group(1)
                contig = match.group(2)
                ## 
                line = fin.readline()
                print("Start", contig)
                pre_dic = {}
                while line:
                    if line.startswith("{corr_start:"):
                        dic = {}
                        ls1 = line.strip("\n").strip(",").lstrip("{").rstrip("}").split(",")
                        for info in ls1:
                            key, val = info.split(":")
                            dic[key] = val
                        line = fin.readline()
                        ls2 = line.strip("\n").strip(",").lstrip("{").rstrip("}").split(",")
                        for info in ls1:
                            key, val = info.split(":")
                            val = val.strip()
                            key = key.strip()
                            dic[key] = val
                        ## 
                        chr_id = dic["chr"]
                        chr_id = int(eval(chr_id))
                        ref_chr = ref_dic[str(chr_id)]
                        if len(pre_dic) > 0:
                            bias = 100
                            # mis_l = int(pre_dic["end"]) - bias
                            # mis_r = int(dic["start"]) + bias
                            # fo.write("{}\t{}\t{}\n".format(ref_chr, str(mis_l), str(mis_r)))
                            # fo.write("{}\t{}\t{}\n".format(ref_chr, dic["start"], dic["end"]))
                            #
                            mis_l = int(pre_dic["end_in_contig"]) - bias
                            mis_r = int(dic["start_in_contig"]) + bias
                            fo.write("{}\t{}\t{}\n".format(contig, str(mis_l), str(mis_r)))
                        pre_dic = dic
                        line = fin.readline()
                    else:
                        break
            elif line.startswith("references_by_id"):
                input_string = line.strip("\n").strip(";")
                regex = r'references_by_id\["([^"]+)"\] = "([^"]+)"'    # references_by_id["0"] = "NC_001133.9"
                match = re.match(regex, input_string)
                print(match.group(1), match.group(2))
                ref_dic[match.group(1)] = match.group(2)
                line = fin.readline()
            else:
                line = fin.readline()

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

# html_file = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/quast/icarus_viewers/alignment_viewer.html"
# bed = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast/quast/icarus_viewers/mis.bed" 
# html_file = sys.argv[1]
# bed = sys.argv[2]

## Test
# html_file = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/quast/icarus_viewers/alignment_viewer.html"
# bed = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/quast/icarus_viewers/mis.bed"

# html_file = sys.argv[1]
# bed = sys.argv[2]
# convert2(html_file, bed)
# convert3(html_file, bed)


# tsv_file = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/quast2/contigs_reports/all_alignments_GCF_000146045.tsv"
# bed = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/yeast2/quast2/contigs_reports/contigs_mis.bed"
if __name__ == "__main__":
    tsv_file = sys.argv[1]
    bed = sys.argv[2]
    bias = int(sys.argv[3])
    convert4(tsv_file, bed, bias)
# python /public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/convert2.py all_alignments*.tsv mis.bed 0