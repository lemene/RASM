import logging
import os
import subprocess
import gzip

prj_dir = os.path.dirname(os.path.abspath(__file__)) 

logger = logging.getLogger()
def enable_logging(log_file="", debug=False, overwrite=True):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)
    logger.addHandler(console_log)
    
    if log_file:
        file_handler = logging.FileHandler(log_file, mode=("a" if  not overwrite else "w"))
        file_handler.setFormatter(log_formatter)
        logger.addHandler(file_handler)

    logger.setLevel(logging.DEBUG)

enable_logging("llog", True)

def safe_make_dir(dir):
    if not os.path.isdir(dir):
        os.makedirs(dir)

def run_command(cmd):
    logger.info("Running: %s" % cmd)
    subprocess.check_call(cmd, shell=True)
    logger.info("Done: %s" % cmd)

def run_mosdepth(bam:str, prefix:str, win_size:int, mapq:int, threads:int=4, ctg:str=None):

    cmd = "mosdepth -Q %d -b %d -t %d -n %s %s" % (mapq, win_size, threads, prefix, bam)
    if ctg != None:
        cmd += " -c %s" % ctg

    run_command(cmd)


def run_samtools_mpileup(bam, output, ctg, ref, mapq):
    cmd = "samtools mpileup -B -q %d -aa -r %s -f %s %s -o %s" % (mapq, ctg, ref, bam, output)
    run_command(cmd)

def split_contig_by_block(ctgs, threads):
    tlen = sum([ctg_len for _, ctg_len in ctgs])
    bsize = max(tlen // (threads *2), 500000)

    for ctg_name, ctg_len in ctgs:
        s = 0
        while s < ctg_len:
            if s + bsize*1.5 < ctg_len:
                e = s + bsize
            else:
                e = ctg_len
            yield ctg_name, ctg_len, s, e
            s = e
    return

class WinIterator(object):
    def __init__(self, length, win_size, stride):
        self.length = length
        self.win_size = win_size
        self.stride = stride
        self.index = 0
        self.num = (self.length - self.win_size + self.stride - 1) // self.stride + 1

    def size(self):
        return self.num
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.index < self.num:
            s = self.stride * self.index
            e = min(self.stride * self.index + self.win_size, self.length)
            self.index += 1
            return s, e
        else:
            raise StopIteration
        

def open_file(fname, mode):
    if fname.endswith(".gz"):
        return gzip.open(fname, mode) 
    else:
        return open(fname, mode)
    

def is_file_newer(files1, files2):
    if type(files1) == str:
        files1_time = os.path.getmtime(files1)
    else:
        files1_time = min([os.path.getmtime(f) for f in files1])
    if type(files2) == str:
        files2_time = os.path.getmtime(files2)
    else:
        files2_time = max([os.path.getmtime(f) for f in files2])

    return files1_time > files1_time

