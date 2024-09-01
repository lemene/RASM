import os
import sys
def get_mummer_script(new_script, ref, query, work_dir):
    raw_script = "/public/home/hpc214712170/shixf/new_code/assembly/Test_code/resolve_err/Pipe/compare_tests/eval_ref/run_mummer.sh"
    VAR_ls = {"ref":ref, "query":query, "work_dir":work_dir}
    with open(raw_script, "r") as fin, open(new_script, "w") as fo:
        for line in fin:
            if line:
                new_line = line
                if line.startswith("#"):
                    fo.write(line)
                    continue
                ls = line.split("=")
                if ls[0] in VAR_ls:
                    new_line = ls[0] + "=" + VAR_ls[ls[0]] + "\n"
                fo.write(new_line)

if __name__ == "__main__":
    ref = os.path.abspath(sys.argv[1])
    query = os.path.abspath(sys.argv[2])
    work_dir = os.path.abspath(sys.argv[3])
    new_script = work_dir + "/run_mummer.sh"
    if not os.path.exists(work_dir): os.makedirs(work_dir)
    get_mummer_script(new_script, ref, query, work_dir)
    