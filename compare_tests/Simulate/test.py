
import pybedtools

bed="/public/home/hpc214712170/Test/mis_detect/simu/cuteSV_simu/data/sim_tra.bed"
bedfile=pybedtools.BedTool(bed)
bedsrtd=bedfile.sort()
# print(bedfile)
# print(bedsrtd)

for j,x in enumerate(bedsrtd):
    print(j, ": ", x)



