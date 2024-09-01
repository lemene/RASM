#!/bin/bash

fa_in=
coverage=40
model=/public/home/hpc214712170/shixf/tools/pbsim3/data/QSHMM-ONT.model
out_dir=`pwd`/fq
prefix=$out_dir/sim_ont_${coverage}x
mkdir -p $out_dir
#Index
samtools faidx $fa_in

## 1、Get ref and asm
#Randomly generate the positions of misassemlies
perl simu_misassembly_posi.pl $fa_in > position.txt

#Remove the redundant positions 
sort -k1V -k2n position.txt | perl move_redundant.pl > position_redun.txt

#Introduce misassemblies by positions. Two FASTA files will be output: one for the reference and one for simulation.
perl simu_misassembly.pl $fa_in

mv *ref.fasta template_ref.fasta
mv *simu.fasta template_simu.fasta

## 2、Get reads
#reads simulation
eval "$(conda shell.bash hook)"
conda activate pbsim3
genome=template_ref.fasta
echo "-------------------------Run pbsim3 for ${genome}-------------------------"
pbsim --strategy wgs --method qshmm --qshmm $model --depth $coverage --genome $genome --prefix $prefix
#The template_simu.fasta is the final assembly with simulated misassemblies, and the template_ref.fasta is the reference.
