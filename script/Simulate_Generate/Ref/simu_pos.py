import random
import os
import sys

base = ['A', 'T', 'C', 'G']
def rand_r(start, end, integer=True):
    if end <= start:
        return None
    if integer:
        return random.randint(start, end)
    else:
        return random.uniform(start, end)

'''def delete_seq(seq, start, length):
    left = seq[:start]
    right = seq[start + length:]
    return left + right

def insert_seq(ins_seq, target_seq, start):
    left = target_seq[:start]
    right = target_seq[start:]
    return left + ins_seq + right

def random_seq(length):
    return ''.join(random.choice(base) for _ in range(length))'''


##########################################################################################

# 获取输入FASTA文件路径
fa = os.path.abspath(sys.argv[1])    # fa input
min_pro = sys.argv[2]
max_pro = sys.argv[3]

# 检查是否存在fai文件
fai = fa + ".fai"
if not os.path.exists(fai):
    raise FileNotFoundError(f'[Error] 找不到fai文件。可以运行"samtools faidx"生成fai文件。')

# parse fa
def parse_fa(fa):
    hash_seq = {}
    id = ''
    with open(fa, 'r') as fa_file:
        for line in fa_file:
            line = line.strip()
            if line.startswith('>'):
                id = line[1:]
            else:
                hash_seq[id] = hash_seq.get(id, '') + line
    return hash_seq


# parse fai
def parse_fai(fai):
    chr_len_dic = {}
    ctg_list = []
    with open(fai, 'r') as fai_file:
        for line in fai_file:
            if line.strip():
                fields = line.split()
                ctg_id = fields[0]
                chr_len_dic[ctg_id] = int(fields[1])
                ctg_list.append(ctg_id)
    return chr_len_dic,ctg_list

chr_len_dic, ctg_list = parse_fai(fai)
# ctg_list = [ctg_id for ctg_id in ctg_list if int(chr_len_dic[ctg_id][0]) >= 500000]

ref_adjust = 0
simu_adjust = 0

# repeats expanded
posi_list = []
for _ in range(50):
    ran_ctg = random.choice(ctg_list)
    ran_position = rand_r(5000, chr_len_dic[ran_ctg] - 6000)
    ins_len = rand_r(1000, 2000)
    posi_list.append([ran_ctg, ran_position, ins_len, "expanded_3"])

# ... Repeat the process for other cases ...

# Inversion
for _ in range(100):
    ran_ctg = random.choice(ctg_list)
    ran_position = rand_r(5000, chr_len_dic[ran_ctg] - 6000)
    ins_len = rand_r(1000, 2000)
    posi_list.append([ran_ctg, ran_position, ins_len, "inv"])

# Insertion
for _ in range(100):
    ran_ctg = random.choice(ctg_list)
    ran_position = rand_r(5000, chr_len_dic[ran_ctg] - 6000)
    ins_len = rand_r(1000, 2000)
    posi_list.append([ran_ctg, ran_position, ins_len, "ins"])

# Deletion
for _ in range(100):
    ran_ctg = random.choice(ctg_list)
    ran_position = rand_r(5000, chr_len_dic[ran_ctg] - 6000)
    ins_len = rand_r(1000, 2000)
    posi_list.append([ran_ctg, ran_position, ins_len, "del"])

for position in posi_list:
    print("\t".join(map(str, position)))


