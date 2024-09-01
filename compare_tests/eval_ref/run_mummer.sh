#!/bin/bash
# mummer -mum -b -c ref.fa query.fa > mummer.mums   # -threads
# mummerplot -x "[0,275287]" -y "[0,265111]" -postscript -p mummer mummer.mums

ref=    # 
query=      # 
work_dir=   # 
threads=40

# 
# mummer -mum -b -c $ref $query -threads $threads > $work_dir/mummer.mums   # -threads
# mummerplot -postscript -p mummer $work_dir/mummer.mums

nucmer -l 100 -c 1000 -d 10 -t 10 --banded -D 5 $ref $query
delta-filter -i 95 -o 95 out.delta > out.best.delta
dnadiff -d out.best.delta
mummerplot out.best.delta --fat -f -png