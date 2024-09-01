#!/bin/bash

SV_type=
raw_dir=data
ref=$raw_dir/reference_genome.fasta
coverage=30
fa_dir=fa_${SV_type}
fq_dir=fq_${SV_type}_${coverage}x
visor_path=/public/home/hpc214712170/singularity/visor_latest.sif

# singularity shell /public/home/hpc214712170/singularity visor_latest.sif
mkdir -p $fq_dir
# SV_ls=(ins del inv tra dup)

# 1、Genome Modification
singularity exec $visor_path \
VISOR HACk -g $ref -b $raw_dir/sim_${SV_type}.bed -o $fa_dir

# 2、Simulated Alignments Generation
# singularity exec $visor_path \
# VISOR LASeR -g $ref -s $fa_dir -b LASeR.bed -o $fq_dir --coverage 30 --threads 40


# WGS simulation
eval "$(conda shell.bash hook)"
conda activate pbsim3

prefix=$fq_dir/sim_ont
model=/public/home/hpc214712170/shixf/tools/pbsim3/data/QSHMM-ONT.model
genome=$fa_dir/h1.fa
pbsim --strategy wgs --method qshmm --qshmm $model --depth $coverage --genome $genome --prefix $prefix

# 

