# ctg_in = "ch1,ch2,ch3"
# ctg_ls = ctg_in.split(",")
# print(ctg_ls)
# from collections import namedtuple

# Region = namedtuple('Region', ["chr_id", "start", "end"])

# region_ls = [Region("chr1", 0, 1000)]
# print(region_ls)

# region_ls1 = [Region]
# print(len(region_ls1), region_ls1)
# for region in region_ls1:
#     print(region)
# region_ls1.append(Region("chr2", 0,1000))
# print(len(region_ls1))
# for region in region_ls1:
#     print(region.chr_id)
# region_ls2 = [Region]

# ls3 = [Region]

# region_ls1.extend(region_ls2)
# print(region_ls1, region_ls2)


## 
import pysam
import time
bam_in = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/resolve_region/solve_del/test/NC_060930.1_hap1.sort.bam"
bam_index = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/resolve_region/solve_del/test/NC_060930.1_hap1.sort.bam.bai"

print(pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index).references)
bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index)
# for 
# start_t = time.time()
# ctg_len = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index).get_reference_length("NC_060930.1")
# cov = bam_reader.count_coverage("NC_060930.1", 0, ctg_len//10)
# end_t = time.time()
# # print(cov)
# print(end_t - start_t)

## time samtools depth -a -J -Q $min_MQ -r $REG ${bam_file} > $out_prefix.depth
reg = "NC_060930.1:20000-20100"
min_MQ = 20
start_t = time.time()
stream = pysam.depth("-a", "-J", "-Q", str(min_MQ), "-r", reg, bam_in)
print(str(stream))
end_t = time.time()
print(end_t - start_t)
ctg_len = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index).get_reference_length("NC_060930.1")
depth_ls = [0] * ctg_len
for line in stream.split("\n"):
    line_ls = line.split()
    if not line_ls:
        continue
    print(line_ls)
    depth_ls[int(line_ls[1]) - 1] = int(line_ls[2])    # 
    
def cal_depth(bam_in, bam_index, ctg):  # cal depth of a contig, return list of depth, length == contig length
    '''time samtools depth -a -J -Q $min_MQ -r $REG ${bam_file} > $out_prefix.depth    
    # del处depth不为0 -J会把deletion算到depth里面去(在deletion处depth不为0)'''
    min_MQ = 20
    ctg_len = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index).get_reference_length(ctg)
    depth_ls = [0] * ctg_len
    depth_stream = pysam.depth("-a", "-J", "-Q", str(min_MQ), "-r", ctg, bam_in)    # depth_stream接收stdout
    for line in depth_stream.split("\n"):
        line_ls = line.split()
        if not line_ls:
            continue
        depth_ls[int(line_ls[1]) - 1] = int(line_ls[2])    # chr_id pos depth   1-based
    return depth_ls

from collections import namedtuple

Region = namedtuple('Region', ["chr_id", "start", "end"])

a = Region("chr1", 1, 10)
print(type(a.end))
ls = [["chr1", 11, 10], ["chr11", 2, 10], ["chr1", 2, 10], ["chr2", 2, 10]] 
ls2 = [Region(reg[0], reg[1], reg[2]) for reg in ls]
ls2.sort(key=lambda res:(res.chr_id, res.start))
print(ls2)