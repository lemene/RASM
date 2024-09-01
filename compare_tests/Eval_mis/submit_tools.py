

import os
import subprocess


def run_cmd_ls(cmd_ls):
    # print("Run: {}".format(" ".join(cmd_ls)))
    subprocess.check_call(" ".join(cmd_ls), shell=True)

# sample_ls = ["sativa_ont", "chm13_ont", "sativa_hifi", "chm13_hifi"]
sample_ls = ["sativa_ont"]
# sample_ls = ["chm13_ont", "sativa_hifi", "chm13_hifi"]
work_dir = "/public/data/biodata/compu_bio/member/shixianfeng/projects/mis/asm"
# tool_ls = ["flye", "wtdbg2", "hifiasm"]
hifi_ls = ["flye", "hifiasm"]
ont_ls = ["flye", "wtdbg2"]
# script_ls = ["run_CRAQ.sh", "run_inspector.sh", "run_GAEP.sh", "run_mypipe.sh"]
# script_ls = ["run_CRAQ.sh", "run_inspector.sh", "run_GAEP.sh"]
script_ls = ["run_mypipe.sh"]
config_ont = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/find_mis_configs/Config-ont.yaml"
config_hifi = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/find_mis_configs/Config-hifi.yaml"

fq = {"chm13_ont": "/public/data/biodata/compu_bio/member/shixianfeng/projects/data/chm13/ont/rel6_30X.fastq.gz",
"chm13_hifi": "/public/data/biodata/compu_bio/member/shixianfeng/projects/data/chm13/hifi/chm13_hifi.fastq.gz",
"sativa_ont": "/public/data/biodata/compu_bio/member/shixianfeng/projects/data/sativa/ont/SRR25241091_40X.fastq",
"sativa_hifi": "/public/data/biodata/compu_bio/member/shixianfeng/projects/data/sativa/hifi/SRR25241090_40X.fastq"}
# print(fq["chm13_hifi"], fq["chm13_ont"])

# exit()
for sample in sample_ls:
    
    if sample.endswith("ont"):
        tool_ls = ont_ls
        data_type = "ont"
        config = config_ont
    else:
        tool_ls = hifi_ls
        data_type = "hifi"
        config = config_hifi
    for tool in tool_ls:
        for script in script_ls:
            raw_script = work_dir + "/" + script
            new_script = work_dir + "/" + sample + "/" + tool + "/" + script    # path/chm13_ont/flye/run_quast.sh
            flag = True
            if os.path.isfile(new_script):
                print("Assembly result {} is existed, check again!!!, Overwrite: yes/no".format(os.path.join(new_script)))
                flag = input()
                if flag == "yes" or flag == "y":
                    flag = True
                else:
                    False
            if flag:
                with open(raw_script, "r") as fin, open(new_script, "w") as fo:
                    for line in fin:
                        if line.startswith("#"):
                            new_line = line
                        elif line.startswith("tool"):
                            new_line = "tool=" + tool + line[5:]
                        elif line.startswith("data_type"):
                            new_line = "data_type=" + data_type + line[10:]
                        elif line.startswith("fastq"):
                            new_line = "fastq=" + fq[sample] + line[6:]
                        elif line.startswith("config"):
                            new_line = "config=" + config + line[7:]
                        else:
                            new_line = line
                        fo.write(new_line)
                print("Ready running for:{}, yes/no | y/n:".format(new_script))
                run_flag = input()
                if flag == "yes" or flag == "y":
                    flag = True
                if flag:
                    print("Start Run: ", new_script)
                    os.chdir(work_dir + "/" + sample + "/" + tool) # 进入其路径
                    cmd_ls = ["/public/home/hpc214712170/submit.sh", script, str(40)]
                    run_cmd_ls(cmd_ls)
                else:
                    print("Skip running: ", new_script)
            else:
                print("Skip: ", new_script)
            
        # script = work_dir + "/" + sample + "/" + tool + "/" + "run_mypipe.sh"
        # run = False
        # print("Ready running:{}, yes/no | y/n:".format(script))
        # flag = input()
        # if flag == "yes" or flag == "y":
        #     os.chdir(work_dir + "/" + sample + "/" + tool)
        #     cmd_ls = ["~/submit.sh", script, "40"]
        #     run_cmd_ls(cmd_ls)
        # else:
        #     print("Cancel!!!")