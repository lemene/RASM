#!/bin/bash

fa_in=

#Index
samtools faidx $fa_in

## 1、Get ref and asm
#Randomly generate the positions of misassemlies
perl simu_misassembly_posi.pl $fa_in > position.txt
echo "----------Simu mis done----------"

#Remove the redundant positions 
sort -k1V -k2n position.txt | perl move_redundant.pl > position_redun.txt
echo "----------Sort mis done----------"

#Introduce misassemblies by positions. Two FASTA files will be output: one for the reference and one for simulation.
perl simu_misassembly.pl $fa_in
echo "----------Simu ref done----------"

mv *ref.fasta template_ref.fasta
mv *simu.fasta template_simu.fasta
#The template_simu.fasta is the final assembly with simulated misassemblies, and the template_ref.fasta is the reference.

