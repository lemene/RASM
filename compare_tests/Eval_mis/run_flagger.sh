#!/bin/bash
work_dir=`pwd`/flagger
data_dir=`pwd`/../simu
data_type=              # hifi/ont/clr
EXPECTED_COVERAGE=      # 40、50
tool_path=/public/home/hpc214712170/singularity/flagger_v0.4.0.sif      # ~/singularity/flagger_v0.2.sif

mkdir -p $work_dir
start=$(date +%s)
# 1、cal depth of  cov
echo "----------------------Start cal depth------------------"
time samtools depth -aa -Q 0 ${data_dir}/${data_type}/aln2simu.sort.bam > ${data_dir}/${data_type}/read_alignment.depth
echo "----------------------Run cal depth dene!!!------------------"

echo "----------------------Start cal cov------------------"
time singularity exec \
--bind $work_dir:$work_dir,$data_dir:$data_dir \
$tool_path \
depth2cov \
-d ${data_dir}/${data_type}/read_alignment.depth \
-f ${data_dir}/Ref/template_ref.fasta.fai \
-o ${work_dir}/read_alignment.cov

echo "----------------------Run cal cov dene!!!------------------"


# 2、cal frequencies of coverages
echo "----------------------Start cal count------------------"
time singularity exec \
--bind $work_dir:$work_dir \
$tool_path \
cov2counts \
-i ${work_dir}/read_alignment.cov \
-o ${work_dir}/read_alignment.counts
echo "----------------------Run cal count dene!!!------------------"

# 3、Coverage Distribution and Fitting The Mixture Model
echo "----------------------Start fittiing model------------------"
time singularity exec \
--bind $work_dir:$work_dir \
$tool_path \
python3 /home/programs/src/fit_gmm.py \
--counts ${work_dir}/read_alignment.counts \
--cov ${EXPECTED_COVERAGE} \
--output ${work_dir}/read_alignment.table
echo "----------------------Run fittiing model dene!!!------------------"

# 4、find_blocks_from_table
echo "----------------------Start find_blocks_from_table------------------"
time singularity exec \
--bind $work_dir:$work_dir \
$tool_path \
find_blocks_from_table \
-c ${work_dir}/read_alignment.cov \
-t ${work_dir}/read_alignment.table \
-p ${work_dir}/out
echo "----------------------Run find_blocks_from_table dene!!!------------------"

end=$(date +%s)
duration=$((end - start))
echo "Total time: $duration seconds."

# out.error.bed out.duplicated.bed out.collapsed.bed
cat $work_dir/out.error.bed $work_dir/out.duplicated.bed $work_dir/out.collapsed.bed > $work_dir/merge.bed
