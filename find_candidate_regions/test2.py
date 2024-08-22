import numpy as np
import subprocess
import pysam
# a = [1,2]
# b =np.array(a)
# print(b)

s= "....,.,.+1C,AGCT-2**"
dic = {".":0, ",":0, "snp":0, }
num = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
sv_len = 0
sv_len_s = ""
for i in s:
    if i == "." or i == ",":
        dic[i] += 1
    elif i == "+":
        sv_len = ""
        
# cmd = ["time samtools mpileup -B -q 20 -aa -d 100 -r NC_000023.11:135000000-136000000 \
#        -f /public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/corrected_ref/reference.fasta \
#        /public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/step1_mapping/aln.sorted.bam > \
#        /public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/pileup/out.test"]
# subprocess.check_call(" ".join(cmd), shell=True)
# subprocess.Popen(" ".join(cmd), shell=True)

# std = pysam.mpileup("-B", "-q", "20", "-aa", "-d", "100", "-r", "NC_000023.11:135000000-136000000", "-f", "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/corrected_ref/reference.fasta", "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/step1_mapping/aln.sorted.bam")
# time samtools mpileup -B -q 20 -aa -d 100 -r NC_000023.11:135000000-136000000 -f /public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/corrected_ref/reference.fasta /public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/step1_mapping/aln.sorted.bam > /public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_ont/my_pipe/pileup/out.test

import os
# print(os.path.join("path","cov","c"))
import heapq
# print("{}".format(" ".join(["1", "2"])))
print(" ".join([]) + "a")