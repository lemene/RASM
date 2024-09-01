#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate syri
tool=flye
data_type=hifi
ref=
asm=data/${tool}_${data_type}.fa
divergence=5 # 5 for hifi, 20 for ONT
work_dir=`pwd`/syri
mkdir -p $work_dir

echo "---------------Start---------------"

# 
nucmer -t48 --maxmatch -c 500 -b 500 -l 100 $ref $asm -p $work_dir/out
echo "---------------Run nucmer done---------------"

# 
delta-filter -m -i 90 -l 100 $work_dir/out.delta > $work_dir/out_m_i90_l100.delta
echo "---------------Run delta-filter done---------------"

show-coords -THrd $work_dir/out_m_i90_l100.delta > $work_dir/out_m_i90_l100.coords
echo "---------------Run show-coords done---------------"
# chroder [-h] [-n NCOUNT] [-o OUT] [-noref] [-F {T,S,B}] coords ref qry
chroder -o $work_dir/out -F S $work_dir/out_m_i90_l100.coords
echo "---------------Run chroder done---------------"
# 1、Using minimap2 for generating alignment. Any other whole genome alignment tool can also be used.
# minimap2 -ax asm${divergence} -t 48 --eqx $ref $asm > $work_dir/out.sam
# echo "---------------Run minimap2 done---------------"
# # 2、
# echo "---------------Start syri---------------"
# syri -c $work_dir/out.sam -r $ref -q $asm -k -F S
# echo "---------------Run syri done---------------"
