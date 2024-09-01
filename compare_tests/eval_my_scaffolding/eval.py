import os
import subprocess


# sample_ls = ["chm13_hifi/my_pipe2", "chm13_ont/my_pipe", "thaliana_ont/my_pipe", "thaliana_hifi/my_pipe", "yeast_ont/my_pipe", "melanogaster_ont/my_pipe"]
# sample_ls = ["yeast_ont_2"]
'''
chm13_hifi/        Mus_musculus_hifi/ reinhardtii_ont/   sativa_ont/        
chm13_ont/         reinhardtii_hifi/  sativa_hifi/       
'''
# sample_ls = ["chm13_hifi/test_params/my_pipe", "chm13_ont/test_params/my_pipe", "Mus_musculus_hifi/test_params/my_pipe"]
sample_ls = ["chm13_hifi", "Mus_musculus_hifi", "reinhardtii_ont", "sativa_ont", "chm13_ont", "reinhardtii_hifi", "sativa_hifi"]
sample_ls = [sample + "/test_params/my_pipe" for sample in sample_ls]
# work_dir="/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests"
# work_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/bench/" # + "/test_params"
work_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/bench"
def get_quast(file1, file2, out_dir):
    try:
        with open(file1, "r") as f1, open(file2, "w") as fo:
            for line in f1:
                if line.startswith("#"):
                    fo.write(line)
                    continue
                new_line = line
                ## 
                if line.startswith("out_dir="):
                    new_line = "out_dir=" + out_dir + "/step5_scaffolding/quast_eval" + "\n"
                if line.startswith("asm="):
                    new_line = "asm=" + out_dir + "/step5_scaffolding/final.fasta" + "\n"
                fo.write(new_line)
        print("Get {} success".format(f2))
    except:
        print("Failed")
for sample in sample_ls:
    print("---------------------------------------------------------------------")
    f1 = work_dir + "/" + sample + "/" + "run_quast.sh"
    f2 = work_dir + "/" + sample + "/" + "step5_scaffolding" + "/" + "run_quast.sh"
    out_dir = work_dir + "/" + sample
    print("out_dir=", out_dir)
    print("f1:", f1)
    print("f2:", f2)
    # print(f1 + "\n" + f2)
    get_quast(f1, f2, out_dir)
    print("Ready running: y/n, yes/no")
    flag = input()
    if flag == "y" or flag == "yes":
        os.chdir(work_dir + "/" + sample + "/" + "step5_scaffolding")
        cmd_ls = ["/public/home/hpc214712170/submit.sh", f2, str(40)]
        subprocess.check_call(" ".join(cmd_ls), shell=True)
        print("Start running for:{}".format(f2))