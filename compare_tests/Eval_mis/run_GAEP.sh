#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate gaep

bam=
asm=
threads=40
out_dir=gaep
data_type=ont
prefix=mis
# config

echo "-----------------------------Start run GAEP-----------------------------"
# eg:gaep pipe -r genome.fasta --lr TGS.fastq -x pb --sr1 NGS_1.fastq --sr2 NGS_2.fastq -t 3 -c config.txt
# gaep bkp -r $asm -i $fastq -t $threads -x $data_type -d $out_dir 
gaep bkp -r $asm -b $bam -t $threads -x $data_type -d $out_dir -o $prefix

