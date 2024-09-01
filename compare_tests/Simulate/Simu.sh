#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate visorenv && echo "visorenv now"
randomregion_r=~/shixf/tools/VISOR/scripts/randomregion.r
work_dir=`pwd`
num=200 # number of intervals
length=2000 # mean intervals length [50000]
standarddev=150 # standard deviation for intervals length [0]
fa_in=Y12.fa
N_reg_bed=$work_dir/N.bed
SV_bed=$work_dir/random.bed
HACk_dir=$work_dir/hack
# var_ls='insertion,deletion,inversion,tandem duplication,inverted tandem duplication'
var_proportion='35:35:10:10:10'    # variants proportions (-r '30:30:40'))

# pipe for easy SV simu
# 1、Creating simple SVs with VISOR
# chrom.dim.tsv: from fai file; exclude_regions.bed: from N region
cut -f1,2 ${fa_in}.fai > $work_dir/chrom.dim.tsv

echo "----------------------Run ${randomregion_r}--------------------------------"
# echo "Rscript $randomregion_r -d $work_dir/chrom.dim.tsv -n $num -l $length -s $standarddev -v $var_ls -r $var_proportion -x $N_reg_bed | sortBed > $work_dir/random.bed"
Rscript $randomregion_r -d $work_dir/chrom.dim.tsv -n $num -l $length -s $standarddev -v 'insertion,deletion,inversion,tandem duplication,inverted tandem duplication' -r $var_proportion -x $N_reg_bed | sortBed > $work_dir/random.bed

# 2、run HACk.py
conda activate mypipe
hack_path=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Simulate/My_simu/HACk.py
echo "----------------------Run HACk.py--------------------------------"
python $hack_path -g $fa_in -b $SV_bed -o $HACk_dir

