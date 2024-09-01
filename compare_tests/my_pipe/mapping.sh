#!/bin/bash

ulimit -n 10000
# reference_fn="/public/data/biodata/compu_bio/member/huangneng/T2T_CHM13/ncbi-genomes-2022-06-13/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
asm_fn=     # 
fastq_fn=
threads=48
prefix=aln
time minimap2 -ax map-ont ${asm_fn} ${fastq_fn} -t $threads | samtools sort -@ $threads > $prefix.sort.bam && samtools index -@ $threads $prefix.sort.bam
# time minimap2 -ax map-ont ${asm_fn} ${fastq_fn} -t $threads | samtools view -bS -@ $threads $prefix.sam | samtools sort -@ $threads -o $prefix.sort.bam && samtools index -@ $threads $prefix.sort.bam