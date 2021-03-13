#!/usr/bin/env python3

import argparse
import edlib
import os
import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
import SDutils
from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-camn")
    parser.add_argument("-iamn")
    parser.add_argument("-o")

    return parser.parse_args()

def get_name(mn, mn_map):
    return mn_map[mn.id][1].split(".")[-1] + "("+ mn.id[3:] + " -- " + mn_map[mn.id][1][7:].split('.')[0] + ":" + str(mn_map[mn.id][0]) + ")"



def get_color(Z, id, CAmn, mn_map):
    c_names = ['black', 'green', 'yellow', 'red', 'firebrick', 'cyan', 'orange', 'chocolate',
               'purple', 'deepskyblue', 'turquoise', 'slategray', 'cornflowerblue', 'blueviolet',
               'violet', 'purple', 'magenta', 'hotpink', 'royalblue', 'azure', 'khaki', 'peru',
               'beige', 'lightgreen', 'thistle', 'oldlace', 'slateblue', 'sienna', 'lime', 'mistyrose',
               'darkturquoise', 'lightcyan', 'olive', 'ghostwhite', 'hotpink', 'tomato', 'darkorange', 'plum',
               'lavender']
    if id < len(CAmn):
        mn = CAmn[id]
        return c_names[int(mn_map[mn.id][1].split('.')[-1])]

    id1 = int(Z[id - len(CAmn)][0])
    id2 = int(Z[id - len(CAmn)][1])

    cl1 = get_color(Z, id1, CAmn, mn_map)
    cl2 = get_color(Z, id2, CAmn, mn_map)

    return cl1 if cl1 == cl2 else "black"
    #
    # if id1 < len(CAmn) and id2 < len(CAmn):
    #     mn1 = CAmn[id1]
    #     mn2 = CAmn[id2]
    #     clid1 = int(mn_map[mn1.id][1].split('.')[-1])
    #     clid2 = int(mn_map[mn2.id][1].split('.')[-1])
    #     if clid1 == clid2:
    #         return c_names[clid1]
    #
    # return "black"


def main():
    args = parse_args()
    CAmn = list(SeqIO.parse(args.camn, "fasta"))
    IAmn = list(SeqIO.parse(args.iamn, "fasta"))
    mn_map = SDutils.map_mn(CAmn, IAmn)
    CAmn.sort(key=lambda x: mn_map[x.id][1].split('.')[-1])
    disis = [[SDutils.get_dist(str(mn1.seq), str(mn2.seq)) for mn2 in CAmn] for mn1 in CAmn]
    df = pd.DataFrame(disis, index=[get_name(mn, mn_map) for mn in CAmn], columns=[get_name(mn, mn_map) for mn in CAmn])
    df.to_csv(args.o)

    Z = linkage(disis, 'single')
    fig = plt.figure(figsize=(25, 10))
    dn = dendrogram(Z, leaf_label_func=lambda x: get_name(CAmn[x], mn_map) if x < len(CAmn) else "",
                    link_color_func=lambda x: get_color(Z, x, CAmn, mn_map))
    plt.show()


if __name__=="__main__":
    main()