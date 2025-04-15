#!/usr/bin/env python3

import os, sys
import argparse
import traceback

def safe_run(cmd):
    print(cmd)
    assert os.system(cmd) == 0

def main(argv):
    parser = argparse.ArgumentParser(description="simulate reads using pbsim3")
    parser.add_argument("ref", type=str)
    parser.add_argument("--depth", type=int, default=40)
    parser.add_argument("--type", type=str, default="clr")
    parser.add_argument("--pbsim", type=str, default="~/tool/pbsim3")
    
#~/tool/pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ~/tool/pbsim3/data/QSHMM-RSII.model --depth 40 --genome ../simu_ref.fasta --prefix hifi_40x --pass-num 10

    try:
        args = parser.parse_args(argv)

        if args.type == "clr":
            cmd = f"{args.pbsim}/src/pbsim \
                --strategy wgs \
                --method qshmm \
                --qshmm {args.pbsim}/data/QSHMM-RSII.model \
                --depth {args.depth} \
                --genome {args.ref}"
            safe_run(cmd)
            pass
        elif args.type == "ont":
            cmd = f"{args.pbsim}/src/pbsim \
                --strategy wgs \
                --method qshmm \
                --qshmm {args.pbsim}/data/QSHMM-ONT-HQ.model \
                --depth {args.depth} \
                --genome {args.ref}"
            safe_run(cmd)
        elif args.type == "hifi":
            cmd = f"{args.pbsim}/src/pbsim \
                --strategy wgs \
                --method qshmm \
                --qshmm {args.pbsim}/data/QSHMM-RSII.model \
                --depth {args.depth} \
                --genome {args.ref} \
                --pass-num 10"
            safe_run(cmd)
        else:
            pass

    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)


if __name__ == "__main__":
    main(sys.argv[1:])