#!/bin/bash

CRAQ_path=/public/home/hpc214712170/shixf/tools/CRAQ/bin
asm=
bam_in=
work_dir=CRAQ
threads=40
mapping_dir=$work_dir/mapping
craq_dir=$work_dir/craq_out
aln_prefix=$mapping_dir/aln
data_type=ont


# mkdir -p $mapping_dir 
#mkdir -p $craq_dir
echo "-----------------------------Start runCRAQ.sh-----------------------------"
# 1、mapping
# time minimap2 -ax map-$data_type -t $threads $asm $fastq > $aln_prefix.sam
# time samtools sort -O BAM -@ $threads $aln_prefix.sam -o $aln_prefix.sorted.bam && samtools index -@ 40 $aln_prefix.sorted.bam
# rm $aln_prefix.sam
# 2、run craq: $ craq  -g assembly.fa -sms SMS_sort.bam
# perl $CRAQ_path/craq -g $asm -sms $aln_prefix.sorted.bam -t $threads -D $work_dir/craq_out
perl $CRAQ_path/craq -g $asm -sms $bam_in -t $threads -D $work_dir/craq_out
