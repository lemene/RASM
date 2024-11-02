#!/bin/bash
# simulate ont reads of the GAEP
#params: --depth 50 --length-min 5000 --length-max 50,000 - --length-mean 20,000

eval "$(conda shell.bash hook)"
conda activate pbsim2

fa_in=../Ref/template_ref.fasta
out_dir=fq
coverage=50
model=/public/home/hpc214712170/shixf/tools/pbsim2/data/P6C4.model # 
prefix=$out_dir/ont_${coverage}x

# WGS simulation
echo "-------------------------Run pbsim2 for ${fa_in}-------------------------"
mkdir -p $out_dir
pbsim --prefix $prefix --hmm_model $model --depth 50 --length-min 5000 --length-max 50000 --length-mean 20000 $fa_in
echo "---------------Run pbsim2 done---------------"
