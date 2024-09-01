#!/bin/bash

source activate flye

genome_size=   # estimated genome size (for example, 5m or 2.6g)
fastq_path=
out_dir=flye
threads=48
time flye --nano-raw $fastq_path -g $genome_size -o $out_dir -t $threads
