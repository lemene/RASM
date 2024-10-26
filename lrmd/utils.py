import logging
import os
import subprocess


logger = logging.getLogger()
def enable_logging(log_file, debug=False, overwrite=True):
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

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


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
        
