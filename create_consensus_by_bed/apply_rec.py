from create_consensus_by_bed import fasta_parser
from create_consensus_by_bed.Utils import Connect_info, Record
# from Utils import Connect_info, Record
# import fasta_parser
def apply_rec_on_ref(rec_ls, ref_seq_in, N_FILL_SIZE):  # N_fill_size
    ''' 提供参考序列一条contig，及对其的所有操作集 '''
    '''per ctg'''
    if len(rec_ls) == 0:
        return [ref_seq_in]
    if len(rec_ls) == 1:
        if rec_ls[0].operation == "keep_ref":   # 摆烂式
            return [ref_seq_in]
        elif rec_ls[0].operation == "replace_with_denovo":  
            seq_out = rec_ls[0].info.rstrip(",").split(",")
            return seq_out

    ## 其他的情形了
    ''' 注意N填充的大小好像会影响到比对结果，最终影响组装结果 '''
    rec_ls.sort(key=lambda rec:int(rec.start))  # 防止有乱序的
    seq_out = []    # 
    offset = 0
    new_seq = ref_seq_in
    max_big_N_fill = 200000     # 对大的区域进行大的N_filling  很有必要（否则会造成一定的影响）
    # big_N_fill_reg_Len = 500000 # patch的长度与估计的长度差值
    for rec in rec_ls:
        start = rec.start + offset
        end = rec.end + offset
        op_ls = rec.operation.split(",")
        info_ls = rec.info.split(",")
        alter_seq = ""

        ## 对info_ls 和 op_ls进行统计，以计算N_fill_size
        N_fill_num = 0
        patch_len = 0
        for op in op_ls:
            if op == "N_fill": N_fill_num += 1
        for info in info_ls:
            patch_len += len(info)  # 
        if N_fill_num > 0:
            differ_len = abs((end - start) - patch_len)
            # if differ_len > big_N_fill_reg_Len:N_fill_size = max_big_N_fill // N_fill_num    # 填充太多是不是也不太好
            # else: N_fill_size = N_FILL_SIZE
            N_fill_size = min((differ_len // (1000 * N_fill_num)) * 1000 + N_FILL_SIZE, max_big_N_fill)
            print("{}:{}-{} N_fill size: {}bp".format(rec.chr_id, rec.start, rec.end, N_fill_size * N_fill_num))

        ## process
        for i, op in enumerate(op_ls):
            if op == "N_fill":
                alter_seq += ("N" * N_fill_size)    # N_fill_size是不定的
            elif op.endswith("patch"):
                alter_seq += info_ls[i]
            else:   # 对于INV等还有待考虑
                # print(rec.rec.operation)
                print(rec.chr_id, rec.start, rec.end, rec.operation, rec.info, "不支持的operation")
                raise(ValueError)   # 待补充？？？
        offset += len(alter_seq) - (rec.end - rec.start)
        new_seq = new_seq[:start] + alter_seq + new_seq[end:]   # seq modify
    seq_out.append(new_seq)
    return seq_out

def apply_rec_on_ref2(chr_id, rec_ls, ref_dic:dict):
    ''' 提供参考序列一条contig，及对其的所有操作集，将处理后的contig写入consensus_fa_dic中 '''
    ''''''
    '''per ctg'''
    ### 
    consensus_fa_dic = {}
    ref_seq_in = ref_dic[chr_id]
    connect_info = Connect_info(chr_id, [], []) # 

    if len(rec_ls) == 0:
        consensus_fa_dic[chr_id] = ref_seq_in
        connect_info.connect_ls.append(chr_id)
        return connect_info, consensus_fa_dic
        # return [ref_seq_in]
    if len(rec_ls) == 1:
        if rec_ls[0].operation == "keep_ref":   # 摆烂式
            consensus_fa_dic[chr_id] = ref_seq_in
            connect_info.connect_ls.append(chr_id)
            return connect_info, consensus_fa_dic
        elif rec_ls[0].operation == "replace_with_denovo":  # 不记录吧，后面通过碎片给记录进去
            return connect_info, consensus_fa_dic
            

    ## 其他的情形了
    ''' 注意N填充的大小好像会影响到比对结果，最终影响组装结果，因此采取先不填充，断开的形式，最后再进行scaffolding '''
    rec_ls.sort(key=lambda rec:int(rec.start))  # 防止有乱序的
    seq_out = []    # 
    # offset = 0
    new_seq = ""
    pre_end = 0    # 
    for rec in rec_ls:  # 重新构写一下，引入对N_fill的处理
        start = rec.start
        end = rec.end
        op_ls = rec.operation.split(",")
        info_ls = rec.info.split(",")
        new_seq += ref_seq_in[pre_end:rec.start]    # 加上上一段rec的末尾到现在rec前面的中间这段序列

        patch_len = 0
        for idx, op in enumerate(op_ls):
            if op.endswith("patch"):patch_len += len(info_ls[idx])  # 求patch length
        ## process
        for i, op in enumerate(op_ls):
            if op == "N_fill":  # 遇到N_fill的位置，直接断开
                if len(new_seq) > 0: 
                    seq_out.append(new_seq)     ## 可能还要记录一下别的信息
                    if info_ls[i] == "INV":
                        connect_info.gap_ls.append(str(end - start) + "_INV")    # 
                        print("INV_region->",chr_id, start, end)
                    else:
                        # print("N_fill_region->",chr_id, start, end, "ref_span:", end - start, "patch_len:", patch_len, "difference:", str(end - start - patch_len) + "bp")
                        print("N_fill_region -> {}:{}-{}, ref_span:{}, patch_len:{}, difference:{}bp".format(chr_id, start, end, end - start, patch_len, end - start - patch_len))
                        connect_info.gap_ls.append(str(end - start - patch_len) + "_GAP")   # 记录下difference or 直接记录end - start距离   ？？
                new_seq = ""    # 更新new_seq
            elif op.endswith("patch"):
                new_seq += info_ls[i]   # 接上这段
                # alter_seq += info_ls[i]
            else:   # 对于INV等还有待考虑
                # print(rec.rec.operation)
                print(rec.chr_id, rec.start, rec.end, rec.operation, rec.info, "不支持的operation")
                # logger.info(rec.chr_id, rec.start, rec.end, rec.operation, rec.info, "不支持的operation")
                raise(ValueError)   # 待补充？？？
        pre_end = rec.end
        

        # offset += len(alter_seq) - (rec.end - rec.start)
        # new_seq = new_seq[:start] + alter_seq + new_seq[end:]   # seq modify
    new_seq += ref_seq_in[pre_end:len(ref_seq_in)]  # 加上最后一个rec后面的序列
    if len(new_seq) > 0: 
        seq_out.append(new_seq)

    ## record to dic
    for i,seq in enumerate(seq_out):
        seq_id = chr_id + "_" + str(i)
        consensus_fa_dic[seq_id] = seq
        connect_info.connect_ls.append(seq_id)  # 按顺序记录seqid的在基因组的连接顺序    
    
    return connect_info, consensus_fa_dic

if __name__ == "__main__":

    chr_id = "NC_000018.10"
    ctg = "NC_000018.10"
    rec_ls = []
    ref_dic = fasta_parser.read_sequence_dict("/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe2/corrected_ref/reference.fasta")
    print("get ref Done")
    # consensus_fa_dic = {}
    consensus_fa_dic = fasta_parser.read_sequence_dict("/public/home/hpc214712170/shixf/projects/ref-guided-assembly/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe/step3_SV_consensus/denovo_fragments.fasta")
    with open("/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe2/step3_SV_consensus/consensus.bed", "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            ctg, start, end, op, info, patch_id = fields[:6]
            # if ctg != chr_id:continue
            start, end = int(start), int(end)
            rec = Record(ctg, start, end)
            rec.add_info(info)
            rec.add_operation(op)
            rec.add_patch_id(patch_id)
            rec_ls.append(rec)
    print("get rec Done")
    connect_info = apply_rec_on_ref2(chr_id, rec_ls, ref_dic, consensus_fa_dic)
    import copy
    fa_out = "/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_18/my_pipe/step3_SV_consensus/test.fa"
    merge_fa_dic = copy.deepcopy(consensus_fa_dic)
    # fasta_parser.write_fasta_dict(consensus_fa_dic, fa_out)
    fasta_parser.write_fasta_dict(merge_fa_dic, fa_out)
    ## 
    # connect_info = Connect_info(chr_id, ["chr1", "chr2", "chr3"], [])
    fo = "/public/home/hpc214712170/Test/tests/chm13_hifi_ctgs/NC_18/test_connect.txt"
    Connect_info.write_connect_info([connect_info], fo)
    pass
