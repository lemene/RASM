# DB_code
These repository includes custom Python scirpts used in publication of DeBreak.


## simulate_sv_genome.py
This script was used to simulate random SVs and embed them into the reference genome.

## sim_find_pos.py
This script was a part of simulation, used to find position of simulated SVs after modification of the genome. 

## extract_sv_vcf.py
This script was used to extract Sv calls located in autosomes and chrX with size >=45bp from VCF file.

## benchmark_callset.py
This script was used to compare two SV callsets (DEL, INS, DUP, and INV). Two SVs are considered as a match if 1) they have same SV type; 2) the breakpoints are located within 1000bp; 3) ratio of SV sizes is less than 2. 

## benchmark_callset_tra.py
This script was used to compare the Translocation calls of two SV callsets. Two TRA calls are considered as a match if both breakpoints are located within 1000bp. 


