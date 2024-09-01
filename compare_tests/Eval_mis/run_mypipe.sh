#!/bin/bash
source activate assembly
script_dir="/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/pipe_main.py"
quast_route="/public/home/hpc214712170/shixf/tools/quast-5.0.2"
work_dir=my_pipe
ref=  # ref for assembly
fq_path=  # fastq in
data_type=	# data type: ont hifi
bam=
genome_size=	# 10m 100m 1G  1000b
ctg_ls=     # ctg_ls kept to perform assembly
threads=40    # 
config=      # config file path
# time python $script_dir -t $threads -a --work-dir $work_dir --ref $ref1 --fastq $fq_path --data-type $data_type -g $genome_size --config $config
/usr/bin/time -v python $script_dir -t $threads --ctg-ls $ctg_ls --work-dir $work_dir --ref $ref --fastq $fq_path --data-type $data_type -g $genome_size --config $config --find_mis --bam $bam

