#!/usr/bin/env python3
 
import sys, os
import traceback
import gzip
import logging
import multiprocessing
import argparse



def safe_run(cmd):
    print(cmd)
    assert os.system(cmd) == 0


def build_bam(ref, rds, bam, threads, params):
    cmd = "minimap2 %s %s -t %d %s -a | samtools sort -@ %d > %s" % (ref, rds, threads, params, threads, bam) 
    print(cmd)
    r = os.system(cmd)
    assert r == 0
    cmd = "samtools index -@ %d %s" % (threads, bam)
    r = os.system(cmd)
    assert r == 0

def md_build_bam(argv):
    '''map reads to reference and generate bamfile'''
    parser = argparse.ArgumentParser(md_build_bam.__doc__)
    parser.add_argument("reference", type=str)
    parser.add_argument("reads", type=str)
    parser.add_argument("bam", type=str)
    parser.add_argument("threads", type=int)
    parser.add_argument("params", type=str)
    try:
        args = parser.parse_args(argv)
        build_bam(args.reference, args.reads, args.bam, args.threads, args.params)

    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)


def md_simu_mis(argv):
    '''simulate misassembly'''
    parser = argparse.ArgumentParser(description=md_simu_mis.__doc__)
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
        get_mis_from_gaep("misassembly.bed", "mis_asm.bed")

    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)


def md_simu_reads(argv):
    '''simulate reads using pbsim3'''
    parser = argparse.ArgumentParser(description=md_simu_reads.__doc__)
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


def get_mis_from_gaep(ifname, ofname):
    with open(ofname, "w") as f:
        for i, line in enumerate(open(ifname)):
            if i % 4 != 3: continue
            its = line.split()
            print(its)
            if its[3] == "del" or its[3] == "ins":
                f.write("%s %s %s\n" % (its[2], its[-2], its[-1]))
            elif its[3] == "inv":
                f.write("%s %s %s\n" % (its[2], its[-4], its[-3]))
            elif its[3] == "collasped":
                f.write("%s %s %s\n" % (its[2], its[-4], its[-3]))
            elif its[3] == "expanded_3":
                f.write("%s %s %s\n" % (its[2], its[-4], its[-3]))
            elif its[3] == "expanded_4":
                f.write("%s %s %s\n" % (its[2], its[-4], its[-3]))
            else:
                print(f"{its[3]}")
                assert not "never come here"


def md_get_mis_from_gaep(argv):
    '''collect misassemblies from gaep output'''
    parser = argparse.ArgumentParser(description=md_simu_mis.__doc__)
    parser.add_argument("ifname", type=str)
    parser.add_argument("ofname", type=str)

    try:
        args = parser.parse_args(argv)
        get_mis_from_gaep(args.ifname, args.ofname)
    except:
        traceback.print_exc()
        print("----------------")
        print(parser.usage)


if __name__ == '__main__':
    if len(sys.argv) > 1:
       locals()[sys.argv[1]](sys.argv[2:])
    else:
       for func in list(locals().keys()):
           if func.startswith("md_"):
               print("%s: %s" % (func, locals()[func].__doc__.split("\n")[0]))
