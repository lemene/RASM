class Feature(object):
    
    @staticmethod
    def calc_win_num(ctg_len, win_size, stride):
        return (ctg_len - win_size + stride - 1) // stride + 1