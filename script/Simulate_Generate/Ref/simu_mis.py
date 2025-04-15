#!/usr/bin/env python3

import sys
import os
import random

def rand_r(start, end, integer=True):
    if end <= start:
        return
    if integer:
        return random.randint(start, end)
    else:
        return random.uniform(start, end)

def random_INV_rev(seq, inv_start, inv_len):
    inv = seq[inv_start:inv_start + inv_len]
    inv = inv[::-1]
    left = seq[:inv_start]
    right = seq[inv_start + inv_len:]
    seq[:] = left + inv + right

def random_INV(seq, inv_start, inv_len):
    inv = seq[inv_start:inv_start + inv_len]
    inv = inv[::-1]
    left = seq[:inv_start]
    right = seq[inv_start + inv_len:]
    seq[:] = left + inv + right

def seq_comp(seq, start, length):
    comp = seq[start:start + length].translate(str.maketrans("ATCGatcg", "TAGCtagc"))
    left = seq[:start]
    right = seq[start + length:]
    seq[:] = left + comp + right

def delete_seq(seq, start, length):
    left = seq[:start]
    right = seq[start + length:]
    seq[:] = left + right

def insert_seq(ins_seq, target_seq, start):
    left = target_seq[:start]
    right = target_seq[start:]
    target_seq[:] = left + ins_seq + right

def random_seq(length):
    bases = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(bases) for _ in range(length))

def random_base_substitute(seq):
    bases = ['A', 'T', 'C', 'G']
    seq_len = len(seq)
    ran_position = random.randint(1, seq_len)
    left = seq[:ran_position - 1]
    right = seq[ran_position:]
    ori_base = seq[ran_position - 1]
    sub_base = random.choice([base for base in bases if base != ori_base])
    seq[:] = left + sub_base + right

def induce_variants(seq, seq_len=None):
    if seq_len is None:
        seq_len = len(seq)
    identity = random.uniform(0.75, 1.0)
    var_threshold = int(seq_len * (1 - identity))

    while var_threshold > 0:
        var_type = random.randint(0, 6)
        if var_type < 5:  # snv
            random_base_substitute(seq)
            var_threshold -= 1
        elif var_type == 5:  # ins
            ran_position = random.randint(1, seq_len)
            ins_len = random.randint(1, 100)
            insert_seq(random_seq(ins_len), seq, ran_position)
            var_threshold -= ins_len
        else:  # del
            ran_position = random.randint(1, seq_len)
            del_len = random.randint(1, 100)
            del_len = min(del_len, seq_len - ran_position)
            delete_seq(seq, ran_position, del_len)
            var_threshold -= del_len

def main():
    fa = sys.argv[1]
    fai = f"{fa}.fai"

    if not os.path.exists(fai):
        sys.exit('[Error] Can\'t find the fai file. You can run "samtools faidx" to generate a fai file.')

    bases = ['A', 'T', 'C', 'G']

    hash_simu = {}
    with open(fa, 'r') as fa_file:
        id = ''
        for line in fa_file:
            line = line.strip()
            if line.startswith('>'):
                id = line[1:]
            else:
                hash_simu[id] = line

    fai_dict = {}
    ctg_list = []
    with open(fai, 'r') as fai_file:
        for line in fai_file:
            if line.strip():
                line_list = line.split()
                fai_dict[line_list[0]] = [int(line_list[1]), int(line_list[2])]
                ctg_list.append(line_list[0])

    ctg_list = [ctg for ctg in ctg_list if fai_dict[ctg][0] >= 500000]

    ref_adjust = 0
    simu_adjust = 0
    prectg = ''
    id_seq = 1
    with open("position_redun.txt", 'r') as pos_file:
        for line in pos_file:
            line = line.strip()
            ctg, posi, length, var_type = line.split()
            if prectg != ctg:
                ref_adjust = 0
                simu_adjust = 0

            posi_ref = int(posi) + ref_adjust
            posi_simu = int(posi) + simu_adjust

            if var_type == 'expanded_3':
                repeat = hash_simu[ctg][posi_ref:posi_ref + int(length)]
                trans_len = rand_r(0, 1000)
                transition = hash_simu[ctg][posi_ref + int(length):posi_ref + int(length) + trans_len]
                repeat_var = induce_variants(repeat, int(length))

                insert_seq(repeat_var, hash_simu[ctg], posi_simu + int(length) + len(transition))
                print('\t'.join([str(id_seq), '1', ctg, var_type, length, str(posi_ref),
                                 str(posi_ref + int(length) + len(transition) + len(repeat_var))]))
                print('\t'.join([str(id_seq), '2', ctg, var_type, length, str(posi_simu),
                                 str(posi_simu + int(length) * 2 + len(transition) * 2 + len(repeat_var) * 2),
                                 str(posi_simu + int(length) + len(transition) + len(repeat_var)),
                                 str(posi_simu + int(length) * 2 + len(transition) * 2 + len(repeat_var))]))

