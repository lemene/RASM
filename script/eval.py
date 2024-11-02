
import os
import sys
from collections import defaultdict
import argparse
import traceback

def load_bed(fname):
    result = []
    for line in open(fname):
        its = line.split()
        if its[0].startswith("#"): continue
            
        chr, start, end = its[:3]
        start = int(start)
        end = int(end)
        result.append((chr, start, end))
    return result


def cal_f1(precision, recall):
    # if precision
    if precision + recall != 0:
        return 2*precision*recall/(precision + recall)
    else:
        return 0
    
def cal_jaccarrd2(comm1, comm2, a, b):
    '''
    将comm1和a放缩为comm2和a_
    '''
    if comm1 == 0 or comm2 == 0: return 0
    a_ = (comm2 / comm1) * a
    return comm2/(comm2 + a_ + b)


def check(truth, test, cov=0.4, bias=0):
    if truth[0] == test[0]:
        assert truth[2] > truth[1] and test[2] > test[1]

        call_reg = [max(test[1]-bias, 0), test[2]+bias]
        inter = min(truth[2], test[2]) - max(truth[1], test[1])
        return inter / (truth[2] - truth[1]) > cov
    else:
        return False


def eval(truth_regions, test_regions, cov,  bias):
    truth_regions.sort()
    test_regions.sort()

    truth_detected = [0] * len(truth_regions)
    test_decteted = [0] * len(test_regions)
    pair = []
    for itruth, truth in enumerate(truth_regions):
        for itest, test in enumerate(test_regions):
            print(truth, test)
            if check(truth, test, cov, bias):
                truth_detected[itruth] += 1
                test_decteted [itest] += 1
                pair.append((itruth, itest))
    
    TP = sum([1 for i in test_decteted if i > 0])
    FP = sum([1 for i in test_decteted if i == 0])
    FN = sum([1 for i in truth_detected if i == 0])

    print(TP, FP, FN)
if __name__ == "__main__":
    parser = argparse.ArgumentParser("")
    parser.add_argument("truth", type=str)
    parser.add_argument("test", type=str)
    parser.add_argument("threshold", type=float)
    parser.add_argument("bias", type=int)


    try:
        args = parser.parse_args(sys.argv[1:])
        truth_regions = load_bed(args.truth)
        test_regions = load_bed(args.test)

        eval(truth_regions, test_regions, args.threshold, args.bias)

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()
