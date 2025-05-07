#!/usr/bin/env python3

import os
import sys
from collections import defaultdict
import argparse
import traceback

class Bed(object):
    def __init__(self, fname, type=""):
        if type == "gaep":
            self.beds = Bed.load_from_gaep(fname)
        else:
            self.beds = Bed.load(fname) 

    def __len__(self):
        return len(self.beds)
    
    def __iter__(self):
        return self.beds.__iter__()
    
    @staticmethod
    def load(fname):
        result = []
        for line in open(fname):
            its = line.split()
            if its[0].startswith("#"): continue
        
            chr, start, end = its[:3]
            result.append((chr, int(start), int(end), line))
        result.sort()
        return result
    def load_from_gaep(fname):
        result = []
        for i, line in enumerate(open(fname)):
            if i % 4 == 3:
                its = line.split()
                if its[3] == 'inv':
                    chr, start, end, = its[2], its[6], its[7]
                else:
                    chr, start, end, = its[2], its[5], its[6]
                result.append((chr, int(start), int(end), line))
                if result[-1][1] >= result[-1][2]:
                    print(result[-1])
                    assert 0
        result.sort()
        return result

class Compare(object):
    def __init__(self, bed_truth, bed_test, threshold:int, bias:int):
        self.truth = bed_truth
        self.test = bed_test
        self.threshold = threshold
        self.bias = bias

        (self.result_truth, self.result_test, self.pairs) = \
            self.compare(bed_truth, bed_test, threshold, bias)

    def compare(self, bed_truth, bed_test, threshold:int, bias:int):
        result_truth = [0] * len(bed_truth)
        result_test = [0] * len(bed_test)
        pairs = []
        for itruth, truth in enumerate(bed_truth):
            for itest, test in enumerate(bed_test):
                if self.check(truth, test, threshold, bias):
                    result_truth[itruth] += 1
                    result_test [itest] += 1
                    pairs.append((itruth, itest))
        return result_truth, result_test, pairs

    @staticmethod
    def check(truth, test, cov=0.4, bias=0):
        if truth[0] == test[0]:
            assert truth[2] >= truth[1] and test[2] >= test[1]
            inter = min(truth[2], test[2]) - max(truth[1], test[1])
            print(truth, test, min(truth[2], test[2]) , max(truth[1], test[1]), inter + bias > 0)
            return inter + bias > 0
        else:
            return False

    def stat(self):
        TP = sum([1 for i in self.result_test if i > 0])
        FP = sum([1 for i in self.result_test if i == 0])
        FN = sum([1 for i in self.result_truth if i == 0])
        
        print(TP, FP, FN)
        precision = TP / (TP+FP)
        recall = TP / (TP + FN)
        print("Precision: %.02f%%" % (precision*100))
        print("Recall: %.02f%%" % (recall*100))
        if precision+recall != 0:
            print("F1: %.02f%%" % (100*(2*precision*recall)/(precision+recall)))
        else:
            print("F1: 0%%")

        c1 = sum([1 for i in self.result_truth if i > 0])
        a1 = sum([1 for i in self.result_truth if i == 0])
        c2 = sum([1 for i in self.result_test if i > 0])
        a2 = sum([1 for i in self.result_test if i == 0])
        print("Jaccard: %.02f%%" % (100*Compare.jaccard(c1, c2, a1, a2)))

    def detail(self):
        print("truth")
        for b, r in zip(self.truth, self.result_truth):
            if r == 0:
                print(b)

        print("test")
        for b, r in zip(self.test, self.result_test):
            if r == 0:
                print(b)
        
    @staticmethod
    def jaccard(c1, c2, a1, a2):
        """a|(c1,c2)|b"""
        if c1 == 0 or c2 == 0: return 0
        a1_ = (c1 / c2) * a1
        return c2/(c2 + a1_ + a2)



    
if __name__ == "__main__":
    parser = argparse.ArgumentParser("")
    parser.add_argument("truth", type=str)
    parser.add_argument("test", type=str)
    parser.add_argument("--truth_type", type=str, default="gaep")
    parser.add_argument("--threshold", type=float, default=0.4)
    parser.add_argument("--bias", type=int, default=500)
    parser.add_argument("--detail", action="store_true", default=False)


    try:
        args = parser.parse_args(sys.argv[1:]) 

        bed_truth = Bed(args.truth, args.truth_type)
        bed_test = Bed(args.test)
        cmp = Compare(bed_truth, bed_test, args.threshold, args.bias)
        cmp.stat()
        if args.detail:
            cmp.detail()

    except:
        traceback.print_exc()
        print("-----------------")
        parser.print_usage()
