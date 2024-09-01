from collections import defaultdict

file = "/public/home/hpc214712170/Test/mis_detect/simu/My_simu/sativa/simu/Ref/mis_simu.bed"
dic = defaultdict(list)
with open(file, "r") as f:
    for line in f:
        ls = line.strip().split("\t")
        ctg, start, end = ls[0], int(ls[1]), int(ls[2])
        # type = ls[4]
        dic[ctg].append(ls)

    for ctg, val_ls in dic.items():
        length = 0
        for ls in val_ls:
            if ls[4] == "ins":
                length += int(ls[6])
            elif ls[4] == "del":
                length += int(ls[6])
            elif ls[4] == "expanded_4" or ls[4] == "expanded_3":
                length += int(ls[2]) - int(ls[1])
            elif ls[4] == "collasped" or ls[4] == "inv":
                length += int(ls[2]) - int(ls[1])
            else:
                raise ValueError
        print(ctg, length)