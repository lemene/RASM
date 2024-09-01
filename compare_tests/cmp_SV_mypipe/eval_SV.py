
'''
思路：
1、提取SV文件
1) vcf过滤
bcftools view -i 'SVTYPE="INS"' cuteSV.vcf.gz -Ou -o cuteSV.ins.vcf     # 按类型过滤
bcftools filter -i 'SVLEN>50 || SVLEN<-50' cuteSV.vcf.gz -Oz -o cuteSV_50.vcf.gz    # 按长度过滤

bcftools view -i 'SVTYPE!="INS" | (SVTYPE=="INS" & SVLEN>50)' input.vcf > output.vcf
bcftools view -i 'SVTYPE=="INV"' cuteSV.vcf > Test.vcf
2) vcf->bed
python vcf2bedpe.py --vcf in.vcf --bedpe out.bed     # vcf -> bedpe
3) 

2、获取mypipe的区域文件，选择consensus.bed作为评估结果
cat consensus.bed | cut -f1-4 > consensus.bed.reg
#chr_id start   end     opertion
NC_000001.11    0       731700  N_fill

3、计算

'''
'''
Commands: 
# bcftools view -i '(SVTYPE=="INS" & SVLEN>50) || (SVTYPE=="DEL" & SVLEN<-50) || (SVTYPE!="INS & SVTYPE!="DEL")' cuteSV.vcf > cuteSV.eval.vcf
# bcftools view -i 'SVTYPE!="INS" & SVTYPE!="DEL" | (SVTYPE=="INS" & SVLEN>50) | (SVTYPE=="DEL" & SVLEN<-50)' cuteSV.vcf > cuteSV.eval.vcf
# python /public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/eval_SV/vcf2bedpe.py --vcf cuteSV.eval.vcf --bedpe cuteSV.eval.bed
'''
'''
选择>100的
'''
'''
# 1、cat cuteSV_100.bed | cut -f1,2,5,11,12 > 
1、cat cuteSV_100.bed | cut -f1,2,5 > out.bed
2、python /public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/eval_SV/eval_SV.py cuteSV_BND.bed cuteSV_BND.cal.bed
3、cat cuteSV_100.cal.bed cuteSV_BND.cal.bed | sort -k 1,1 -k 2n,2 | less > merge.bed
4、bedtools merge -d 100 -i merge.bed > merge.out.bed
'''

# chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstarnd1\tstrand2\tsvtype\tnumber_of_support_read

import sys
import os

def process_BND(bed_in, bed_out):
    if bed_in == bed_out: raise ValueError
    with open(bed_in, "r") as infile, open(bed_out, "w") as outfile:
        for line in infile:
            if line.startswith("#"):continue
            # print(line)
            parts = line.strip().split("\t")
            chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, svtype, number_of_support_read = parts

            # 提取第一个区间并扩展50bp
            new_start1 = max(int(start1) - 50, 0)
            new_end1 = int(end1) + 50

            # 提取第二个区间并扩展50bp
            new_start2 = max(int(start2) - 50, 0)
            new_end2 = int(end2) + 50

            # 创建新的行并写入输出文件
            new_line = f"{chrom1}\t{new_start1}\t{new_end1}\n{chrom2}\t{new_start2}\t{new_end2}\n"
            outfile.write(new_line)


fin = os.path.abspath(sys.argv[1])
fout = os.path.abspath(sys.argv[2])
process_BND(fin, fout)

