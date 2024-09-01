''' split asm to eval'''
import os
import subprocess
import sys
sys.path.append("/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe")
from create_consensus_by_bed.Utils import make_dir, Record, run_cmd_ls
from create_consensus_by_bed import fasta_parser
from compare_tests.collect_data_of_quast import collect_specify
import pysam
from multiprocessing import Pool
import shutil

def get_quast_script(quast_file, new_quast_file, asm, Reference, out_dir, threads, large, quast_params):
    with open(quast_file, "r") as fin, open(new_quast_file, "w") as fo:
        for line in fin:
            if fin:
                new_line = line
                if line.startswith("#"):
                    fo.write(line)
                    continue
                if line.startswith("asm="):
                    new_line = "asm=" + asm + line[len("asm="):]
                if line.startswith("Reference="):
                    new_line = "Reference=" + Reference + line[len("Reference="):]
                if line.startswith("out_dir="):
                    new_line = "out_dir=" + out_dir + line[len("out_dir="):]
                if line.startswith("threads="):
                    new_line = "threads=" + str(threads) + line[len("threads="):]
                if line.startswith("large="):
                    new_line = "large=" + large + line[len("large="):]
                if line.startswith("params="):
                    new_line = "params=" + "\"" + str(quast_params) + "\"" + line[len("params="):]
                fo.write(new_line)
    
def run_submit(script, args):
    '''Usage: sbatch [OPTIONS(0)...] [ : [OPTIONS(N)...]] script(0) [args(0)...]'''
    cmd_ls = ["sbatch -D `pwd` -A pi_zy -J", script, "--partition cpuQ -q cpuq", "--cpus-per-task=48", "-o", script+".log", script, args]
    run_cmd_ls(cmd_ls)
    pass

def quast_asm(asm_ls, quast_dir, Reference, quast_params_ls):
    threads = 48
    make_dir(quast_dir)
    quast_script_path = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Eval/run_quast.sh" 
    for asm in asm_ls:
        name = os.path.basename(asm).rstrip(".fa")
        quast_out_dir = quast_dir + "/" + name
        new_quast_script_path = quast_dir + "/run_quast_" + name + ".sh"
        print("Get script:", new_quast_script_path)
        quast_params = " ".join(([str(arg) for arg in quast_params_ls]))
        get_quast_script(quast_script_path, new_quast_script_path, asm, Reference, quast_out_dir, threads, large, quast_params)
        args_ls = [asm, Reference, quast_out_dir, threads, large, quast_params]
        print("args: ", " ".join([str(arg) for arg in args_ls]))
        run_submit(new_quast_script_path, "")

def cut_seq(cut_info_ls, seq_in:str, ctg):
    ctg_patch_dic = {}
    ctg_consensus_dic = {}
    id = 0
    pre_end = 0
    for cut_info in cut_info_ls:
        start, end, patch_id = cut_info
        patch_seq = seq_in[start:end]
        consensus_seq = seq_in[pre_end:start]
        print("{}, seq: {}-{}, {}-{}".format(ctg, pre_end, start, start, end), end="; ")
        if len(consensus_seq) > 0:
            ctg_consensus_dic[ctg + "_" + str(id)] = consensus_seq
        if len(patch_seq) > 0:
            # ctg_patch_dic[patch_id] = patch_seq
            ctg_patch_dic[ctg + "_frage_" + str(id)] = patch_seq
        pre_end = end
        id += 1
    consensus_seq = seq_in[pre_end:]
    if len(consensus_seq) > 0:
        ctg_consensus_dic[ctg + "_" + str(id)] = consensus_seq
        print("{}, seq: {}-{}".format(ctg, pre_end, len(seq_in)))
    return ctg_patch_dic, ctg_consensus_dic

def _split(bam_in, ref_fa_dic:dict):
    '''组装几部分：
    1、fragemant
    2、reads_patch (skip!!!)
    3、asm_patch
    4、consensus
    '''
    # frage_dic = {}
    consensus_dic = {}
    patch_dic = {}
    all_info = []
    bam_reader = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_in + ".bai")
    for ctg in bam_reader.references:   # 遍历ctg，将ctg在对应位置切开
        ctglen = bam_reader.get_reference_length(ctg)
        ctg_seq = ref_fa_dic[ctg]
        cut_info = []
        for read in bam_reader.fetch(ctg, 0, ctglen):
            if read.is_secondary or read.is_unmapped: continue # suppl ?
            if read.is_supplementary: 
                # print("{} is suppl, alignlen:{}".format(read.query_name, read.query_alignment_length))
                continue
            ref_start, ref_end = read.reference_start, read.reference_end
            if ref_start < 100: ref_start = 0 
            if ref_end > ctglen - 100: ref_end = ctglen
            patch_id = read.query_name
            cut_info.append([ref_start, ref_end, patch_id])
        ctg_patch_dic, ctg_consensus_dic = cut_seq(cut_info, ctg_seq, ctg)
        patch_dic.update(ctg_patch_dic), consensus_dic.update(ctg_consensus_dic)
        all_info.append([ctg, ctg_patch_dic, ctg_consensus_dic, cut_info])
    for info in all_info:
        ctg, ctg_patch_dic, ctg_consensus_dic, cut_info = info
        if len(cut_info) > 0:
            print("Cut {} -> {}, {}".format(ctg, [new_ctg for new_ctg in ctg_patch_dic.keys()], [new_ctg for new_ctg in ctg_consensus_dic.keys()]))
            print(ctg, "cut_info:", cut_info)
    return patch_dic, consensus_dic


def eval_mypipe(out_dir, consensus_bed, merge_fn, polish_fn, large, Reference, quast_params_ls, run_q):
    make_dir(out_dir)
    # 1、get split asm from rec and merge.fasta/ro polish
    rec_ls = Record.read_record(consensus_bed)
    replace_with_denovo = []
    asm_patch = []
    N_fill = []
    asm_patch_dic1 = {}
    for rec in rec_ls:
        if rec.operation == "reads_patch": continue
        elif rec.operation == "replace_with_denovo": replace_with_denovo.append(rec)
        elif rec.operation == "asm_patch": 
            asm_patch.append(rec)
            asm_patch_dic1[rec.patch_id] = rec.info
        else:
            op_ls = rec.get_op_ls()
            info_ls = rec.get_info_ls()
            patch_ls = rec.get_patch_ls()
            if rec.operation == "N_fill":
                N_fill.append(rec)
                continue
            elif "asm_patch" in op_ls:
                asm_patch.append(rec)
                for idx, op in enumerate(op_ls):
                    if op == "asm_patch":
                        if patch_ls[idx] in asm_patch_dic1: print("{}:{}-{} use dup".format(rec.chr_id, rec.start, rec.end))
                        asm_patch_dic1[patch_ls[idx]] = info_ls[idx]
            else: raise ValueError
    
    polish_dic = fasta_parser.read_sequence_dict(polish_fn)
    frage_dic = {}
    consensus_dic = {}
    for ctg, seq in polish_dic.items():
        if ctg.endswith("fragemant"): frage_dic[ctg] = seq
        else: consensus_dic[ctg] = seq
    
    consensus_fn1 = out_dir + "/polished.consensus.fa"
    asm_patch_fn1 = out_dir + "/raw.asm_patch.fa"
    frage_fn = out_dir + "/polished.frage.fa"
    fasta_parser.write_fasta_dict(consensus_dic, consensus_fn1)
    fasta_parser.write_fasta_dict(asm_patch_dic1, asm_patch_fn1)
    fasta_parser.write_fasta_dict(frage_dic, frage_fn)
    # mapping patch_seq
    bam_out = out_dir + "/" + "patch_to_consensus.bam"
    minimap2_cmd_ls = ["minimap2", "-ax", "asm20", "-t20", consensus_fn1, asm_patch_fn1, \
        "|", "samtools", "sort", "-O", "BAM", "-@20", "-o", bam_out, \
        "&&", "samtools index -@20", bam_out] 
    run_cmd_ls(minimap2_cmd_ls)
    # _split
    for rec in asm_patch:
        print("{}:{}-{}, {}, {}".format(rec.chr_id, rec.start, rec.end, rec.operation, rec.patch_id))
    asm_patch_fn2 = out_dir + "/polished.asm_patch.fa"
    consensus_fn2 = out_dir + "/polished.expatch.consensus.fa"
    patch_dic, polish_consensus_dic= _split(bam_out, consensus_dic)
    fasta_parser.write_fasta_dict(patch_dic, asm_patch_fn2)
    fasta_parser.write_fasta_dict(polish_consensus_dic, consensus_fn2)
    ## 2、run eval
    threads = 40
    asm_ls = [consensus_fn1, consensus_fn2, asm_patch_fn1, asm_patch_fn2, frage_fn]
    quast_dir = out_dir + "/quast"
    quast_script_path = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Eval/run_quast.sh" 
    if run_q:
        quast_asm(asm_ls, quast_dir, Reference, quast_params_ls)
    print("my_asm_ls: ", asm_ls)
    return asm_ls

def run_eval_mypipe(work_dir, large, Reference, quast_params_ls):
    eval_dir = work_dir + "/split_eval"
    consensus_bed = work_dir + "/step3_SV_consensus/consensus.bed"
    merge_fn = work_dir + "/step3_SV_consensus/merge.fasta"
    polish_fn = work_dir + "/step4_polish/racon.fasta" 
    asm_ls = eval_mypipe(eval_dir, consensus_bed, merge_fn, polish_fn, large, Reference, quast_params_ls, True)
    # out_dir, consensus_bed, merge_fn, polish_fn, large, Reference, quast_params_ls, run_q
    return asm_ls
def eval_other(out_dir, frage_fn, asm_fn, Reference, quast_params_ls):
    make_dir(out_dir)
    ## 
    # mapping patch_seq
    bam_out = out_dir + "/" + "frage_to_asm.bam"
    minimap2_cmd_ls = ["minimap2", "-ax", "asm20", "-t30", asm_fn, frage_fn, \
        "|", "samtools", "sort", "-O", "BAM", "-@30", "-o", bam_out, \
        "&&", "samtools index -@30", bam_out] 
    run_cmd_ls(minimap2_cmd_ls)
    ## 
    asm_dic = fasta_parser.read_sequence_dict(asm_fn)
    tool_frage_dic, tool_consensus_dic = _split(bam_out, asm_dic)
    tool_frage_fn = out_dir + "/" + "tool_frage.fa"
    tool_consensus_fn = out_dir + "/" + "tool_consensus.fa"
    fasta_parser.write_fasta_dict(tool_frage_dic, tool_frage_fn)
    fasta_parser.write_fasta_dict(tool_consensus_dic, tool_consensus_fn)
    ## 
    threads = 48
    quast_dir = out_dir + "/quast"
    make_dir(quast_dir)
    quast_script_path = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/Eval/run_quast.sh" 
    asm_ls = [tool_frage_fn, tool_consensus_fn]
    print("tool_asm_ls: ", asm_ls)
    return asm_ls

def collect_report(quast_dir, name_ls):
    print("\nQuast report:")
    dir_Ls = os.listdir(quast_dir)
    for dir in dir_Ls:
        report_file_path = quast_dir + "/" + dir + "/report.txt"
        if not os.path.isfile(report_file_path): continue
        print(report_file_path)
        collect_specify(report_file_path, name_ls)

if __name__ == "__main__":

    my_work_dir = "/public/home/hpc214712170/Test/tests/thaliana_ont/my_pipe3"
    large = "yes"
    Reference = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/thaliana/REF/ref2/GWHBDNP00000000.genome.fasta.gz"
    quast_params_ls = ["--min-identity", "90.0"]
    my_quast_dir = my_work_dir + "/split_eval/quast"
    # my_asm_ls = run_eval_mypipe(my_work_dir, large, Reference, quast_params_ls)
    my_asm_ls = ['/public/home/hpc214712170/Test/tests/thaliana_ont/my_pipe3/split_eval/polished.consensus.fa', '/public/home/hpc214712170/Test/tests/thaliana_ont/my_pipe3/split_eval/polished.expatch.consensus.fa', '/public/home/hpc214712170/Test/tests/thaliana_ont/my_pipe3/split_eval/raw.asm_patch.fa', '/public/home/hpc214712170/Test/tests/thaliana_ont/my_pipe3/split_eval/polished.asm_patch.fa', '/public/home/hpc214712170/Test/tests/thaliana_ont/my_pipe3/split_eval/polished.frage.fa']
    # quast_asm(my_asm_ls, my_work_dir + "/split_eval/quast", Reference, quast_params_ls)

    tool_out_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/flye/split_eval"
    frage_fn = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/my_pipe3/split_eval/polished.asm_patch.fa"
    asm_fn = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/flye/assembly.fasta"
    Reference = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/data/thaliana/REF/ref2/GWHBDNP00000000.genome.fasta.gz"
    quast_params_ls = ["--min-identity", "90.0"]
    tool_quast_dir = tool_out_dir + "/quast"
    # tool_asm_ls = eval_other(tool_out_dir, frage_fn, asm_fn, Reference, quast_params_ls)
    tool_asm_ls = ['/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/flye/split_eval/tool_frage.fa', '/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/thaliana_ont/flye/split_eval/tool_consensus.fa']
    # quast_asm(tool_asm_ls, tool_out_dir + "/quast", Reference, quast_params_ls)
    ## 

    name_ls = ["NG50", "NGA50", "Genome fraction (%)", "indels per 100 kbp", "mismatches per 100 kbp", "QV", "misassemblies", "local misassemblies"]
    name_ls = ["Genome fraction (%)", "misassemblies", "local misassemblies"]
    collect_report(my_quast_dir, name_ls)
    collect_report(tool_quast_dir, name_ls)
