

import os
import subprocess


def run_cmd_ls(cmd_ls):
    # print("Run: {}".format(" ".join(cmd_ls)))
    subprocess.check_call(" ".join(cmd_ls), shell=True)
sample_ls = ["sativa_ont", "chm13_ont", "sativa_hifi", "chm13_hifi"]
work_dir = "/public/home/hpc214712170/Test/mis_detect/asm"
tool_ls = ["flye", "wtdbg2", "hifiasm"]
for sample in sample_ls:
    for tool in tool_ls:
        script = work_dir + "/" + sample + "/" + tool + "/" + "run_mypipe.sh"
        run = False
        print("Ready running:{}, yes/no | y/n:".format(script))
        flag = input()
        if flag == "yes" or flag == "y":
            os.chdir(work_dir + "/" + sample + "/" + tool)
            cmd_ls = ["~/submit.sh", script, "40"]
            run_cmd_ls(cmd_ls)
        else:
            print("Cancel!!!")