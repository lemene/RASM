import os
import sys
quast_script = "run_quast.sh"
def convert_size(size:str):
    if size.endswith("G") or size.endswith("g"):
        return 1000000*float(size[:-1])
    elif size.endswith("M") or size.endswith("m"):
        return 1000*float(size[:-1])
    elif size.endswith("B") or size.endswith("b"):
        return float(size[:-1])
    else:
        raise ValueError("Error format of genome size!!!")
def get_quast_script(quast_file_in, out_dir, ref1, ref2, genome_size):
    quast_file_out = os.path.join(out_dir, quast_script)
    with open(quast_file_in, "r") as fin, open(quast_file_out, "w") as fo:
        for line in fin:
            if fin:
                new_line = line
                if line.startswith("#"):
                    fo.write(line)
                    continue
                if line.startswith("asm="):
                    new_line = "asm=" + ref1 + line[len("asm="):]
                if line.startswith("Reference="):
                    new_line = "Reference=" + ref2 + line[len("Reference="):]
                if line.startswith("out_dir="):
                    new_line = "out_dir=" + out_dir + "/quast_eval" + line[len("out_dir="):]
                if line.startswith("large="):
                    genome_size_num = convert_size(genome_size)
                    if genome_size_num > 100000:
                        new_line = "large=" + "yes" + line[len("large="):]
                    else:
                        new_line = "large=" + "no" + line[len("large="):]
                fo.write(new_line)


if __name__ == "__main__":
    quast_script_in = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/修改模板脚本并提交任务/eval_ref/run_quast.sh"
    out_dir = os.path.abspath(sys.argv[1])
    ref1 = os.path.abspath(sys.argv[2])
    ref2 = os.path.abspath(sys.argv[3])
    g_size = sys.argv[4]
    print("\nout_dir="+out_dir)
    print("ref1="+ref1)
    print("ref2="+ref2)
    print("genome_size="+g_size)
    try:
        get_quast_script(quast_script_in, out_dir, ref1, ref2, g_size)
    except:
        print("Error")
    pass