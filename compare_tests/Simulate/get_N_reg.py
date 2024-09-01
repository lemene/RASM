# 提取序列中的N区域

import sys
# from Bio import SeqIO
fa_in = sys.argv[1]
bed_out = sys.argv[2]
def parse_fa(fa):
    seq_dic = {}
    id = ''
    with open(fa, 'r') as fa_file:
        for line in fa_file:
            line = line.strip()
            if line.startswith('>'):
                id = line[1:]
            else:
                seq_dic[id] = seq_dic.get(id, '') + line
    return seq_dic

def extract_N_regions(fa_in, bed_out):
    fo = open(bed_out, "w")
    seq_dic = parse_fa(fa_in)
    min_length = 50
    for ctg, seq in seq_dic.items():
        n_regions = []
        current_region = None
        for idx, base in enumerate(seq):
            if base == "N":
                if current_region is None:
                    current_region = [idx, idx]
                else:
                    current_region[1] = idx
            elif current_region is not None:
                if current_region[1] - current_region[0] + 1 >= min_length:
                    n_regions.append(current_region)
                current_region = None
        if current_region is not None and current_region[1] - current_region[0] + 1 >= min_length:
            n_regions.append(current_region)
        # write bed
        print(ctg + ":" ,n_regions)
        for reg in n_regions:
            fo.write("{}\t{}\t{}\n".format(ctg, reg[0], reg[1]))
    fo.close()

extract_N_regions(fa_in, bed_out)

