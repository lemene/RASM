#!/bin/bash
source activate assembly
## 
ctg_ls=
ref=
bam_in=
threads=
config_file=
out_dir=    # pileup

script_dir=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/find_candidate_regions/find_from_pileup.py
time python $script_dir $ctg_ls $ref $bam_in $threads $config_file $out_dir