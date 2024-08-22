import sys
import numpy as np
# import scipy.stats as st
file = sys.argv[1]
reverse_flag = sys.argv[2]
reverse_flag = int(reverse_flag)
ls = []
win_size = 1000
str_ls= []
None_val_num = 0
val_num = 0
with open(file, "r")as f:
    print("Start")
    for line in f:
        if line.startswith("#"): 
            print(line)
            continue
        fields = line.strip()
        val_num += 1
        if not fields:
            None_val_num += 1
            continue
        ls.append(float(line.strip()))
    print("None_val porotion: ", None_val_num / val_num)
    print("mean: ", np.mean(ls))
    print("median: ", np.median(ls))
    step = 0.001
    quantile_ls = [i * step for i in range(int(1 / step))]
    quantile_num_ls = np.quantile(ls, quantile_ls)
    if reverse_flag:
        quantile_num_ls = np.flip(quantile_num_ls)
    for i in range(len(quantile_ls)):
        print("{:.3f}:{:.7f}".format(quantile_ls[i], quantile_num_ls[i]))

def stats():
    return
def run():
    pass