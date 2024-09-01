#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate ins
inspector_path=/public/home/hpc214712170/shixf/tools/Inspector

# # Evaluate assembly with raw reads
# inspector.py -c contig.fa -r rawreads.1.fastq rawreads.2.fastq -o inspector_out/ --datatype clr 
# # Evaluate assembly with hifi reads
# inspector.py -c contig.fa -r ccsreads.1.fastq ccsreads.2.fastq -o inspector_out/ --datatype hifi

fastq=
asm=
out_dir=inspector
data_type=

echo "-----------------------------Start run inspector.py-----------------------------"
$inspector_path/inspector.py -c $asm -r $fastq -o $out_dir --datatype $data_type --min_contig_length_assemblyerror 100

