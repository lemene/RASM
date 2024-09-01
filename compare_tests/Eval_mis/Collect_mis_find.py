

import os
import subprocess
import time


def collect_stats(file, val_ls):
    dic = dict.fromkeys(val_ls)
    if os.path.isfile(file):
        with open(file, "r") as f:
            for line in f.readlines():
                fields = line.strip().split(":")
                # print(fields)
                dic[fields[0]] = fields[1].strip()
                # dic
            for val in val_ls:
                print(dic[val], end="|")
            print("\n")
    else:
        print(None)
        return

def run_cmd_ls(cmd_ls):
    t0 = time.time()
    # print("Run: {}".format(" ".join(cmd_ls)))
    subprocess.check_call(" ".join(cmd_ls), shell=True)
    # print("Run: {} finished, cost {}s".format(" ".join(cmd_ls), time.time() - t0))
if __name__ == "__main__":
    sample_ls = ["sativa_ont", "chm13_ont", "sativa_hifi", "chm13_hifi"]
    # work_dir = "/public/home/hpc214712170/Test/mis_detect/asm"
    work_dir = "/public/data/biodata/compu_bio/member/shixianfeng/projects/mis/asm"
    tool_ls = ["flye", "wtdbg2", "hifiasm"]
    val_ls = ["Contigs", "Length", "Longest contig", "N50", "QV"]
    print(val_ls)
    for sample in sample_ls:
        print("---------------------------------{}--------------------------------".format(sample))
        for tool in tool_ls:
            print("********************************{}*******************************".format(tool))
            target_dir = work_dir + "/" + sample + "/" + tool + "/my_pipe"
            stats_file = target_dir + "/step2_candidate_regions/simple.ststs"
            collect_stats(stats_file, val_ls)
            ## 
            print("------Mis num-------")
            mis_bed = target_dir + "/step2_candidate_regions/filtered2/merge/merge.bed"
            # medge_ = target_dir + "/step2_candidate_regions/filtered/merge.bed"
            # medge_ = target_dir + "/step2_candidate_regions/filtered2/merge.bed"
            # if os.path.isfile(medge_):
            #     cmd_ls = ["cat", medge_, "| wc -l"]
            #     run_cmd_ls(cmd_ls)
            # else:
            #     # mis_call_dir = target_dir + "/step2_candidate_regions/filtered/*"
            #     # mis_call_dir = target_dir + "/step2_candidate_regions/filtered2/*"
            #     # cmd_ls = ["cat", mis_call_dir, "| wc -l"]
            #     run_cmd_ls(cmd_ls)
            if os.path.isfile(mis_bed):
                cmd_ls = ["wc -l", mis_bed]
                run_cmd_ls(cmd_ls)
            else:
                print(None)
