#At first, we recommend to remove the contigs less than 500k in template.fasta and split the template.fasta from Ns.

#Index
samtools faidx template.fasta

#Randomly generate the positions of misassemlies
perl simu_misassembly_posi.pl template.fasta > position.txt

#Remove the redundant positions 
sort -k1V -k2n position.txt | perl move_redundant.pl > position_redun.txt

#Introduce misassemblies by positions. Two FASTA files will be output: one for the reference and one for simulation.
perl simu_misassembly.pl template.fasta

mv *ref.fasta template.ref.fasta
mv *simu.fasta template_simu.fasta

#PacBio reads simulation
eval "$(conda shell.bash hook)"
conda activate visorenv

echo "Run pbsim2"
pbsim template.ref.fasta --prefix simu_pb --depth 50 --length-min 5000 --length-max 50000 --hmm_model ~/shixf/tools/pbsim2/data/P6C4.model --length-mean 20000

#The template_simu.fasta is the final assembly with simulated misassemblies, and the template.ref.fasta is the reference.
