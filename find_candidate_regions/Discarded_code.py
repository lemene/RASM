# def find_by_cov(bam_in, ctg, bed_out, MIN_DP, MIN_DP_REG):
#     bam_index = bam_in + ".bai"
#     ## cal_depth
#     print("\n------------------cal depth of {}------------------".format(ctg))
#     dp_ls = cal_depth(bam_in, ctg)   # cal depth of this contig
#     ctg_len = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index).get_reference_length(ctg)
#     AVG_DP = sum(dp_ls) / ctg_len
#     print("AVG_DP:", AVG_DP)    # AVG_DP: 26

#     # MIN_DP = 3
#     # MIN_DP_REG = 100   # 500
#     win_size = 100
#     win_num = ctg_len // win_size + 1 if ctg_len % win_size != 0 else ctg_len // win_size
#     # print(ctg_len % win_size != 0)
#     avg_dp_ls = [0] * win_num
#     low_dep_reg = []
#     for i in range(win_num):
#         if i < win_num - 1:
#             avg_dp_ls[i] = cal_sum(dp_ls, i * win_size, (i + 1) * win_size) / win_size
#             if avg_dp_ls[i] < MIN_DP:
#                 low_dep_reg.append(Region(ctg, i*win_size, (i+1)*win_size))    # append low_dep
#         else:  # ==win_num-1
#             avg_dp_ls[i] = cal_sum(dp_ls, i * win_size, ctg_len) / (ctg_len - i * win_size)
#             if avg_dp_ls[i] < MIN_DP:
#                 low_dep_reg.append(Region(ctg, i*win_size, ctg_len))    # append low_dep
#         # print("{}:{}-{}  dp:{}".format(ctg, i * win_size, (i + 1) * win_size, avg_dp_ls[i]))
#         ## 
#         # if avg_dp_ls[i] < MIN_DP:
#         #     low_dep_reg.append(Region(ctg, i*win_size, (i+1)*win_size))    # append low_dep
#     # print(low_dep_reg)
#     low_dep_reg = cluster_regions2(low_dep_reg)
#     # low_dep_reg_merge.extend(low_dep_reg)
#     ##
#     # bed_out = os.path.join(out_dir, ctg+"cov.bed")
#     with open(bed_out, "w") as fout:
#         for reg in low_dep_reg:
#             if reg.end - reg.start >= MIN_DP_REG:   # 记录大于500bp的
#                 fout.write("{}\t{}\t{}\t{}bp-low_dep\n".format(reg.chr_id, reg.start, reg.end, reg.end-reg.start))
#                 # fout.write("{}\t{}\t{}\t{}\n".format(ctg, reg[0], reg[1], "low_dep"))
#             print("{}\t{}\t{}".format(reg.chr_id, reg.start, reg.end))

# def find_by_cov2(bam_in, ctg, bed_out, MIN_DP, MIN_DP_REG):
#     bam_index = bam_in + ".bai"
#     ## cal_depth
#     print("\n------------------cal depth of {}------------------".format(ctg))
#     dp_ls = cal_depth(bam_in, ctg)   # cal depth of this contig
#     ctg_len = pysam.AlignmentFile(bam_in, "rb", index_filename=bam_index).get_reference_length(ctg)
#     AVG_DP = sum(dp_ls) / ctg_len
#     print("AVG_DP:", AVG_DP)    # AVG_DP: 26

#     # MIN_DP = 3
#     # MIN_DP_REG = 100   # 500
#     win_size = 100
#     win_num = ctg_len // win_size + 1 if ctg_len % win_size != 0 else ctg_len // win_size
#     # print(ctg_len % win_size != 0)
#     avg_dp_ls = [0] * win_num
#     low_dep_reg = []
#     for i in range(win_num):
#         if i < win_num - 1:
#             avg_dp_ls[i] = cal_sum(dp_ls, i * win_size, (i + 1) * win_size) / win_size
#         else:  # ==win_num-1
#             avg_dp_ls[i] = cal_sum(dp_ls, i * win_size, ctg_len) / (ctg_len - i * win_size)
#         # print("{}:{}-{}  dp:{}".format(ctg, i * win_size, (i + 1) * win_size, avg_dp_ls[i]))
#         ## 
#         if avg_dp_ls[i] < MIN_DP:
#             low_dep_reg.append(Region(ctg, i*win_size, (i+1)*win_size))    # append low_dep
#     # print(low_dep_reg)
#     low_dep_reg = cluster_regions2(low_dep_reg)
#     # low_dep_reg_merge.extend(low_dep_reg)
#     ##
#     # bed_out = os.path.join(out_dir, ctg+"cov.bed")
#     with open(bed_out, "w") as fout:
#         for reg in low_dep_reg:
#             if reg.end - reg.start >= MIN_DP_REG:   # 记录大于500bp的
#                 fout.write("{}\t{}\t{}\t{}bp-low_dep\n".format(reg.chr_id, reg.start, reg.end, reg.end-reg.start))
#                 # fout.write("{}\t{}\t{}\t{}\n".format(ctg, reg[0], reg[1], "low_dep"))
#             print("{}\t{}\t{}".format(reg.chr_id, reg.start, reg.end))
#     return