#!/usr/bin/env python3

import argparse
import edlib
import os
import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
import SDutils
from matplotlib import pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-sum")
    parser.add_argument("-mon")
    parser.add_argument("-o")

    return parser.parse_args()

def handle_cen(cenid, args, df):
    cenid_set = {"'" + cenid + "0'", "'" + cenid + "1'", "'" + cenid + "2'", "'" + cenid + "3'"}
    print(cenid_set)
    print(set(df[0][-6]))
    df_cen = [x for x in df if len(cenid_set & set(x[-6])) > 0]
    mn_id_sum = []
    for i in range(len(df_cen)):
        sum = 0
        for j in range(len(df_cen[i][-6])):
            if df_cen[i][-6][j] in cenid_set:
                sum += int(df_cen[i][-5][j])
        mn_id_sum.append([df_cen[i][0], sum])
    mn_id_sum.sort(key=lambda x: -x[1])

    mons = list(SeqIO.parse(args.mon, "fasta"))
    dmon = {mn.id: mn for mn in mons}
    print(mn_id_sum)

    mn_id_sum[0].append(100)
    mn_id_sum[0].append("-")
    for i in range(1, len(mn_id_sum)):
        cur_sep = 100
        b_mn = "-"
        for j in range(i):
            nm1 = "mn_" + str(mn_id_sum[i][0])
            nm2 = "mn_" + str(mn_id_sum[j][0])
            #print(nm1, nm2)
            sep = SDutils.seq_identity(str(dmon[nm1].seq), str(dmon[nm2].seq))
            if sep < cur_sep:
                cur_sep = sep
                b_mn = nm2
        mn_id_sum[i].append(cur_sep)
        mn_id_sum[i].append(b_mn)

    df = pd.DataFrame(mn_id_sum, columns=["MnId", "Cnt", "Sep", "ClosestMn"])
    df.to_csv(os.path.join(args.o, cenid + ".csv"))

def main():
    args = parse_args()
    df = pd.read_csv(args.sum).values.tolist()
    print(df[0])
    for i in range(len(df)):
        df[i][-6] = df[i][-6].strip('][').split(', ')
        df[i][-5] = df[i][-5].strip('][').split(', ')

    for i in range(1, 23):
        handle_cen("cen"+str(i)+"_", args, df)
    handle_cen("cenX_", args, df)


if __name__ == "__main__":
    main()