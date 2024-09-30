from find_candidate_regions.find_mis_pipe import *

class Detector:
    def __init__(self):
        pass

    

    def run_find_pipe(self, ref, bam, ctgs, out_dir, threads, config):

        dp_win_size = config["dp_params"]["dp_win_size"]
        block_size = config["dp_params"]["block_size"]
        min_MQ = config['dp_params']['min_MQ']
        dp_info_dic = get_dp_info_parallel(bam, threads, out_dir, dp_win_size, block_size, min_MQ)

        ## ----------------1、find candidate----------------
        # dirs
        candidate_dir = os.path.join(out_dir, "candidate")
        if not os.path.isdir(candidate_dir):os.mkdir(candidate_dir)
        info_bed_dir = os.path.join(out_dir, "info")
        if not os.path.isdir(info_bed_dir):os.mkdir(info_bed_dir)
        filtered2_dir = os.path.join(out_dir, "filtered2")
        if not os.path.isdir(filtered2_dir):os.mkdir(filtered2_dir)
        pileup_dir = os.path.join(out_dir, "pileup")
        if not os.path.isdir(pileup_dir):os.mkdir(pileup_dir)

        # run find
        pool = Pool(processes=threads)
        for ctg in ctgs:
            pool.apply_async(find_candidate, args=(ctg[0], 0, ctg[1], ref, bam, out_dir, dp_info_dic[ctg[0]], config))
        
        pool.close() 
        pool.join()

