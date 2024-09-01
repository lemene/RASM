#!/bin/bash
WORKDIR=/public/data/biodata/compu_bio/member/shixianfeng/projects/mis/asm
ex_size=40000    # 200000
cluster_len=5000
data_type_ls=("hifi" "ont")
# sample_ls=("chm13_hifi","chm13_ont","sativa_hifi","sativa_ont")
sample_ls=("chm13" "sativa")
# data_type_ls=("hifi")
# sample_ls=("sativa")
tool_ls=("gaep" "craq" "inspector" "my_pipe")

# echo ${sample_ls[@]}
for sample in "${sample_ls[@]}"
do
    echo "####################### Start $sample #######################"
    for data_type in "${data_type_ls[@]}"
    do
        echo $data_type
        echo "---------------------sample:${sample}_${data_type}---------------------"
        # get asm tools
        if [ "$data_type" = "ont" ]; then
            asm_tool_ls=("flye" "wtdbg2")
        else
            asm_tool_ls=("flye" "hifiasm")
        fi
        # 
        for asm_tool in "${asm_tool_ls[@]}" 
        do
            echo "-----------$asm_tool----------"
            # sample_dir=$WORKDIR/$sample/$data_type
            sample_dir=$WORKDIR/${sample}_${data_type}/$asm_tool
            quast_dir=${sample_dir}/quast-lg
            # 一、For QUAST
            # 1、Generate QUAST mis bed, run convert2.py
            script1_dir=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/convert2.py
            bias=0
            # echo "step1"
            python $script1_dir ${quast_dir}/contigs_reports/all_alignments*.tsv ${quast_dir}/contigs_reports/mis.bed $bias > ${quast_dir}/contigs_reports/mis.bed.log

            # 2、Run mis_clusters.py for QUAST mis bed
            script2_dir=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/mis_cluetsr.py
            # echo "step2"
            python $script2_dir ${quast_dir}/contigs_reports/mis.bed ${quast_dir}/contigs_reports/mis_cl.bed $cluster_len

            # 3、Run ex.py
            script3_dir=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/Get_ex.py
            cat ${quast_dir}/contigs_reports/all_alignments*.tsv | grep CONTIG | cut -f2,3 | sort -n -k 2 > ${quast_dir}/contigs_reports/ctg.lengths
            # echo "step3"
            python $script3_dir ${quast_dir}/contigs_reports/ctg.lengths ${quast_dir}/contigs_reports/ex.txt $ex_size > ${quast_dir}/contigs_reports/ex.txt.log

            # 二、For tool
            echo "Start process tool"
            eval_script=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/Eval_mis_find.py
            threshold=0
            eval_bias=0
            
            for tool in "${tool_ls[@]}"
            do
                echo ">>>${tool}"
                if [ "$tool" = "my_pipe" ]; then
                    mis_bed=${sample_dir}/my_pipe/step2_candidate_regions/filtered2/merge/merge.bed
                elif [ "$tool" = "craq" ]; then
                    cat ${sample_dir}/craq_out/runAQI_out/locER_out/out_final.CRE.bed ${sample_dir}/craq_out/runAQI_out/strER_out/out_final.CSE.bed > ${sample_dir}/craq_out/mis.bed
                    mis_bed=${sample_dir}/craq_out/mis.bed
                elif [ "$tool" = "gaep" ]; then
                    mis_bed=${sample_dir}/gaep/mis_breakpoints.txt
                elif [ "$tool" = "inspector" ]; then
                    mis_bed=${sample_dir}/inspector/structural_error.bed
                else
                    echo "Error tool name"
                fi

                if [ ! -e $mis_bed ]; then
                    echo "-------Failed:${mis_bed}-------"
                else
                    echo "Eval ${mis_bed}"
                    python $eval_script ${quast_dir}/contigs_reports/mis_cl.bed $mis_bed $threshold $eval_bias ${quast_dir}/contigs_reports/ex.txt > ${mis_bed}.log
                    jacc=$(cat ${mis_bed}.eval | grep jacc)
                    echo "${jacc} ${tool}"
                fi
            done
        done
    done
done

exit

echo "---------------------${sample} ${data_type}---------------------"
quast_dir=${sample_dir}/quast-lg
# 一、For QUAST
# 1、Generate QUAST mis bed, run convert2.py
script1_dir=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/convert2.py
bias=0
python $script1_dir ${quast_dir}/contigs_reports/all_alignments*.tsv ${quast_dir}/contigs_reports/mis.bed $bias

# 2、Run mis_clusters.py for QUAST mis bed
script2_dir=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/mis_cluetsr.py
python $script2_dir ${quast_dir}/contigs_reports/mis.bed ${quast_dir}/contigs_reports/mis_cl.bed $cluster_len

# 3、Run ex.py
script3_dir=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/Get_ex.py
cat ${quast_dir}/contigs_reports/all_alignments*.tsv | grep CONTIG | cut -f2,3 | sort -n -k 2 > ${quast_dir}/contigs_reports/ctg.lengths
python $script3_dir ${quast_dir}/contigs_reports/ctg.lengths ${quast_dir}/contigs_reports/ex.txt $ex_size

# 二、For tool
eval_script=/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/Eval__mis_SV/Eval_mis_find.py
threshold=0
eval_bias=0

for tool in tool_ls
do
    echo "-------${tool}-------"
    if [ "$tool"="my_pipe" ]; then
    mis_bed=${sample_dir}/my_pipe/step2_candidate_regions/filtered2/merge/merge.bed
    elif [ "$tool"="craq" ]; then
        cat ${sample_dir}/craq_out/runAQI_out/locER_out/out_final.CRE.bed ${sample_dir}/craq_out/runAQI_out/strER_out/out_final.CSE.bed > ${sample_dir}/craq_out/mis.bed
        mis_bed=${sample_dir}/craq_out/mis.bed
    elif [ "$tool"="gaep" ]; then
        mis_bed=${sample_dir}/gaep/mis_breakpoints.txt
    elif [ "$tool"="inspector" ]; then
        mis_bed=${sample_dir}/inspector/structural_error.bed
    else
        echo "Error tool name"
    fi

    if [ ! -e $mis_bed ]; then
        echo "-------${tool} is none, failed!!!-------"
    else
        python $eval_script ${quast_dir}/contigs_reports/mis_cl.bed $mis_bed $threshold $eval_bias ${quast_dir}/contigs_reports/ex.txt
    fi
done







