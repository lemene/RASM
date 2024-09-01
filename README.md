# LRMD
A genome assembly evaluating pipeline.

## Introduction

LRMD is a reference-free misassembly detection tool that accurately and comprehensively identifies misassembly areas in assembly.

## Installation
```shell
# 
git clone https://github.com/sxfss/LRMD.git
# install the conda envirament
cd LRMD/env
conda env create -n LRMD -f LRMD_env.yml
# Ensure that these tools are installed in your system and added to the environment.
minimap2(we use version 2.24)
samtools(we use version 1.17)
mosdepth(https://github.com/brentp/mosdepth)
```

## Quick Start
***Before running `LRMD` please do not forget to activate the conda environment.***

### Step 1: Create a config file and 
You can directly use Config-ont.yaml and Config-hifi.yaml under the LRMD/find_mis_configs directory as configuration files. 
You can also copy a template Config.yaml file from LRMD\find_mis_configs and then modify the parameters. A complete config file is shown as follows.
```shell
step1:
  minimap2_params_ls: []    # 
step2:
  win_size: 400
  apply_dp: true
  dp_params:
    dp_upper_bound: 1.75   # dp * upper_bpund -> upper_dp
    dp_lower_bound: 0.15    # dp * lower_bpund -> lower_dp
    dp_win_size: 100  # win_size of tool(mosdepth) to calculate depth
    block_size: 10000000   # 
    cluster_dis: 1000   # 
    min_MQ: 0  # 0 for asm   
  apply_clip: true
  clip_params:
    min_clip_len: 500 # 
    min_clip_portion: 0.35  # 
    min_clip_num: 5   # 
    clip_win_size: 1000 #
    min_MQ: 0  # 
  apply_pileup: true
  pileup_params:
    reg_size: 1000000
    win_size: 400   # 
    step_size: 200  # 
    min_correct_portion: 0.9    # 
    max_differ_portion: 0.2     # 
    max_disagree_portion: 0.1  # 
    cluster_dis: 1000
    min_MQ: 0  # 
  filter2:
    bound_span: 10000 #  
    check_win_size: 8000 # 
    min_span_num: 2   # 
    min_supp_portion: 0.4  # 
    min_MQ: 1    # 
    MIN_ALIGN_LENGTH: 10000
    MIN_ALIGN_RATE: 0.95
    ins_threshold: 100
    del_threshold: 100
    min_clip_len: 500
```

### Step 2: Mapping reads
```shell
# Mapping reads to the assembly result.
minimap2 -ax map-ont ${asm} ${fq} -t 40 | samtools sort -@ 40 > aln.sort.bam
samtools index -@ 40 aln.sort.bam
```Some other statistical indicators.

### Step 3: Running LRMD
```shell
python $PATH/LRMD/pipe_main.py -t 40 -a --work-dir LRMD --ref $ref --data-type ont --config Config.yaml --find_mis --bam aln.sort.bam
```

## Output
### Misassembly results
LRMD/step2_candidate_regions/filtered2/merge/merge.bed is the final misassembly results.
### Some statistical indicators
LRMD/step2_candidate_regions/simple.ststs is some other statistical indicators.
