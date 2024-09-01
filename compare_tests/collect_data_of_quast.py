# import csv
import os
import subprocess
import sys
import math
'''
script for collect evaluate value in quast

'''
def cal_QV(mismatches, indels):
    '''QV = 10 * log10(100000/(mismatches + indels))'''
    QV = 10 * math.log(100000.0/(float(mismatches) + float(indels)), 10)
    return QV

def collect(file_ls, file_out, tool_name, sample_ls, name_ls, split):
    # name_ls = ["NG50", "NGA50", "Genome fraction (%)", "indels per 100 kbp", "mismatches per 100 kbp", "misassemblies", "local misassemblies"]
    fo = open(file_out, "w")
    # fo.write("#" + "|".join(name_ls) + "\n")
    # print("#" + "|".join(name_ls))
    fo.write("#{}\n".format(tool_name))
    fo.write("#" + "|".join(name_ls) + "\n")
    for idx, file in enumerate(file_ls):
        dic = {}
        if os.path.isfile(file):
            with open(file, "r")as fin:
                fin.readline()
                fin.readline()
                fin.readline()
                
                for line in fin:
                    line = line.strip()
                    # print(line)
                    if line.startswith("#"):
                        name = line[2:29]
                    else:
                        name = line[:29]
                    info = line[29:].strip()
                    name = name.strip()
                    # print(name)
                    # print(info)
                    if name in name_ls:
                        # print(name)
                        dic[name] = info
        if dic.get("indels per 100 kbp", None) != None:
            dic["QV"] = cal_QV(dic["mismatches per 100 kbp"], dic["indels per 100 kbp"])
        fo.write("\n# {}\n".format(sample_ls[idx]))
        # for name in name_ls:
        #     fo.write("{}\t{}\n".format(name, dic.get(name, None)))
        for name in name_ls:
            fo.write("{}{}".format(dic.get(name, None), split))
    fo.close()

def run(sample_ls, fin_dir, fout_dir, tool, name_ls):
    '''批量统计，并输出统计结果至文件中'''
    fin_ls = []
    for sample in sample_ls:
        quast_file1 = os.path.join(fin_dir, sample, tool, "quast_eval", "report.txt")
        quast_file2 = os.path.join(fin_dir, sample, tool, "quast_eval_large", "report.txt")
        if os.path.isfile(quast_file1):
            fin_ls.append(quast_file1)
        else:
            fin_ls.append(quast_file2)
    for file in fin_ls:
        print(file)
    fout = fout_dir + "/" + tool + "_quast.stats"
    collect(fin_ls, fout, tool, sample_ls, name_ls)

def collect_specify(file, name_ls, split):
    '''进行单个的统计'''
    print("#" + "|".join(name_ls))
    dic = {}
    if os.path.isfile(file):
        with open(file, "r")as fin:
            fin.readline()
            fin.readline()
            fin.readline()
            
            for line in fin:
                line = line.strip()
                # print(line)
                if line.startswith("#"):
                    name = line[2:29]
                else:
                    name = line[:29]
                info = line[29:].strip()
                name = name.strip()
                # print(name)
                # print(info)
                if name in name_ls:
                    # print(name)
                    dic[name] = info
    if dic.get("indels per 100 kbp", None) != None:
        dic["QV"] = cal_QV(dic["mismatches per 100 kbp"], dic["indels per 100 kbp"])
    # for name in name_ls:
    #     print("{}".format(dic.get(name, None)))
    for name in name_ls:
        print("{}".format(dic.get(name, None)), end=split)
    print("\n")
# def merge(tool_ls, dir):
#     cmd_ls = ["paste"]
#     tool_path = [out_dir + "/" + tool + "_quast.stats" for tool in tool_ls]
#     ls1 = [">", dir + "/merge.stats"]
#     cmd_ls.extend(tool_path)
#     cmd_ls.extend(ls1)
#     # print(tool_ls, tool_path)
#     cmd = " ".join(cmd_ls)
#     print(cmd)
#     # print(" ".join(cmd_ls))
#     subprocess.check_call(cmd, shell=True)
#     pass


if __name__ == "__main__":
    # print(cal_QV(2.23, 1.99))
    # exit(0)
    split = "|"
    name_ls = ["Total length", "contigs", "NG50", "NGA50", "Genome fraction (%)", "indels per 100 kbp", "mismatches per 100 kbp", "QV", "misassemblies", "local misassemblies"]
    # # sample_ls = ["yeast_ont", "melanogaster_ont", "chm13_hifi", "chm13_ont", "thaliana_hifi", "thaliana_ont"]
    # sample_ls = ["reinhardtii_ont", "sativa_ont", "chm13_ont", "reinhardtii_hifi", "sativa_hifi", "Mus_musculus_hifi", "chm13_hifi"]
    # work_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests"
    # fin_ls = []
    # # tool = "flye"
    # # tool = "wtdbg2"
    # # tool = "hifiasm"
    # # # tool = "polish"
    # # tool = "my_pipe"
    # # tool = "shasta"
    # # tool = "raven"
    # # tool = "canu"
    # # tool = "smartdenovo"
    # tool = "my_pipe"
    # # work_dir = "/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/bench"
    # for sample in sample_ls:
    #     quast_file1 = os.path.join(work_dir, sample, tool, "quast_eval", "report.txt")
    #     quast_file2 = os.path.join(work_dir, sample, tool, "quast_eval_large", "report.txt")
    #     # quast_file1 = os.path.join(work_dir, sample, "test_params", tool, "quast_eval", "report.txt") #for my_pipe
    #     # quast_file2 = os.path.join(work_dir, sample, "test_params", tool, "quast_eval_large", "report.txt")
    #     if os.path.isfile(quast_file1):
    #         fin_ls.append(quast_file1)
    #     else:
    #         fin_ls.append(quast_file2)
    # # fin_ls = ["/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/melanogaster_ont/my_pipe/quast_eval_large/report.txt"]
    # # fout = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/quast.stats"
    # fout = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/eval_asm" + "/" + tool + "_quast.stats"
    # print("\t".join(sample_ls))
    # for file in fin_ls:
    #     print("Collect:", file)
    # print("Stats Out:", fout)
    # collect(fin_ls, fout, tool, sample_ls, name_ls, split)
    # exit(0)

    ##-----------------------------------------------------------------------------------------------
    # out_dir = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/quast_stats"
    # out_dir = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/eval_asm"
    # tool_ls = ["wtdbg2", "flye", "hifiasm", "shasta", "my_pipe", "polish"]
    # # tool_ls = ["shasta"]
    # name_ls = ["Total length", "contigs", "NG50", "NGA50", "Genome fraction (%)", "indels per 100 kbp", "mismatches per 100 kbp", "QV", "misassemblies", "local misassemblies"]
    # print("\n" + "#" + "|".join(name_ls) + "\n")
    # for tool in tool_ls:
    #     # run(sample_ls, work_dir, out_dir, tool, name_ls)  # 批量计算
    #     pass
    # out_ls = [out_dir + "/" + tool + "_quast.stats" for tool in tool_ls]
    # # merge(tool_ls, out_dir)   # paste f1 f2 > merge
    

    ## collect_specify
    file = os.path.abspath(sys.argv[1])
    collect_specify(file, name_ls, split)


    # collect_specify("/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi/my_pipe2/quast_eval_large/report.txt")
    pass
    