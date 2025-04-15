#!/usr/bin/env python3

import os, sys
import argparse
import traceback

def safe_run(cmd):
    print(cmd)
    assert os.system(cmd) == 0

def main(argv):
    parser = argparse.ArgumentParser(description="simulate misassembly")
    parser.add_argument("ref", type=str)
    parser.add_argument("--gaep", type=str)

    try:
        args = parser.parse_args(argv)

        cmd = f"samtools faidx {args.ref}"
        safe_run(cmd)

        cmd = f"perl {args.gaep}simu_misassembly_posi.pl {args.ref} > position.txt"
        safe_run(cmd)

        cmd = f"sort -k1V -k2n position.txt | perl {args.gaep}move_redundant.pl > position_redun.txt"
        safe_run(cmd)

        cmd = f"perl {args.gaep}simu_misassembly.pl {args.ref} > misassembly.bed"
        safe_run(cmd)

        cmd = f"mv {args.ref}*ref.fasta simu_ref.fasta"
        safe_run(cmd)

        cmd = f"mv {args.ref}*simu.fasta simu_asm.fasta"
        safe_run(cmd)

        # 
        with open("mis_asm.bed", "w") as f:
            for i, line in enumerate(open("misassembly.bed")):
                if i % 4 == 3:
                    its = line.split()
                    f.write("%s %s %s\n" % (its[2], its[-2], its[-1]))

    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)


if __name__ == "__main__":
    main(sys.argv[1:])