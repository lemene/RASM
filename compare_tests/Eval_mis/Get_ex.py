
import sys
def read_fai(fai_filename):
    chromosome_lengths = {}
    try:
        with open(fai_filename, 'r') as fai_file:
            for line in fai_file:
                parts = line.strip().split('\t')  # Split each line into parts using tab as a delimiter
                if len(parts) >= 2:
                    chrom_name = parts[0]  # The first part is the chromosome name
                    chrom_length = int(parts[1])  # The second part is the chromosome length
                    chromosome_lengths[chrom_name] = chrom_length
    except IOError:
        print(f"Error: File {fai_filename} does not exist or cannot be read.")
    return chromosome_lengths

if __name__ == "__main__":

    fai = sys.argv[1]
    ex_f = sys.argv[2]
    ctg_len_dic = read_fai(fai)
    min_ctg = int(sys.argv[3])   # 200000
    skip_ls = []
    ex_num = 0
    num = 0
    for ctg, ctg_len in ctg_len_dic.items():
        if ctg_len_dic[ctg] < min_ctg:
            print("Skip:{}".format(ctg))
            skip_ls.append(ctg)
            ex_num += 1    
        num += 1
        
    with open(ex_f, "w") as f:
        f.write("# all num: {}\n".format(num))
        f.write("# ex num: {}\n".format(ex_num))
        for ctg in skip_ls:
            f.write("{}\n".format(ctg))

# /public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/Get_ex.py

# cat all_alignments*.tsv | grep CONTIG | cut -f2,3 | sort -n -k 2 > ctg.lengths
# python /public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/Get_ex.py ctg.lengths ex.txt 200000
# 