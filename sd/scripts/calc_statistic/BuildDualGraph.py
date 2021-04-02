#!/usr/bin/env python3

import argparse
import edlib
import os
import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
import csv
import math
from matplotlib import pyplot as plt
from subprocess import check_call
import networkx as nx
from networkx.algorithms import bipartite

import TriplesMatrix
from TriplesMatrix import calc_mn_order_stat
import SDutils

def parse_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-sdtsv")
    parser.add_argument("-sep", default="-")
    parser.add_argument("-mon", default="-")
    parser.add_argument("-o")

    return parser.parse_args()

def save_vert(fw, mn1, mn2, cenName, cnt_mon1, cnt_mon2):
    fw.write("graph " + cenName + " {\n")
    fw.write(" "* 4 + "rankdir=LR;\n")

    for mn, cnt_mon in [(mn1, cnt_mon1), (mn2, cnt_mon2)]:
        for vert in mn:
            cnt_mn = 0
            lg = 0.01
            if vert in cnt_mon:
                cnt_mn = cnt_mon[vert]
                lg = math.log(cnt_mn)
            clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa"]
            #print(int(lg))
            fw.write(" " * 4 + '"' + vert + "\" [style=filled fillcolor=\"" + clr[int(lg)] + "\" label=\"" + vert + "[" + str(cnt_mn) + "]\"];\n")

    fw.write(" "*4 + "{rank = same; " + "; ".join(mn1) + ";}\n")
    fw.write(" "*4 + "{rank = same; " + "; ".join(mn2) + ";}\n")


def save_vert_mn(fw, mn1, cenName, cnt_mon1):
    fw.write("digraph " + cenName + " {\n")
    #fw.write(" "* 4 + "rankdir=LR;\n")
    fw.write(" " * 4 + "ratio = 1.0;\n")

    for mn, cnt_mon in [(mn1, cnt_mon1)]:
        for vert in mn:
            cnt_mn = 0
            lg = 0.01
            if vert in cnt_mon:
                cnt_mn = cnt_mon[vert]
                lg = 0
                if cnt_mn > 0:
                    lg = math.log(cnt_mn)
            clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa"]
            #print(int(lg))
            #vertn = "-".join([x[3:] for x in list(vert)])
            vert = "-".join(list(vert))
            fw.write(" " * 4 + '"' + vert + "\" [style=filled fillcolor=\"" + clr[int(lg)] + "\" label=\"" + vert + "[" + str(cnt_mn) + "]\"];\n")


def save_edges(fw, mn1, mn2, db_cnt):
    for vt1 in sorted(mn1):
        fw.write(" " * 4 + "\"" + vt1 + "\" -- \"" + vt1 + "_" + "\"\n")
        for vt2 in sorted(mn2):
            if (vt1, vt2) not in db_cnt:
                continue
            scr = db_cnt[(vt1, vt2)]

            thr_wg = [100000000, 1000, 500, 100, 1]
            wgs = [7, 5, 3, 1, 0]
            wg = 3
            while (scr > thr_wg[wg]):
                wg -= 1

            if scr > 100:
                fw.write(" " * 4 + "\"" + vt1 + "\" -- \"" + vt2 + "\"")
                fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=" + str(wgs[wg]) + "];\n")
            else:
                pass
                #fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=" + str(wgs[wg]) + " constraint=false];\n")

    fw.write("}\n")


def save_edges_mn(fw, mn1, kcnt, matching):
    for vt1 in sorted(mn1):
        for vt2 in sorted(mn1):
            #print(vt1, vt2)
            #print(list(vt1)[1:], list(vt2)[:-1])
            if list(vt1)[1:] != list(vt2)[:-1]:
                continue
            scr = 0
            if (*vt1, vt2[-1]) in kcnt:
                scr = kcnt[(*vt1, vt2[-1])]

            thr_wg = [100000000, 1000, 500, 100, 1]
            wgs = [7, 5, 3, 1, 0]
            wg = 3
            while (scr > thr_wg[wg]):
                wg -= 1

            if matching.get(vt1, "") == (*vt2, "_"):
                vrt1 = "-".join(list(vt1))
                vrt2 = "-".join(list(vt2))
                fw.write(" " * 4 + "\"" + vrt1 + "\" -> \"" + vrt2 + "\"")
                fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=" + str(wgs[wg]) + " color=blue];\n")
            elif scr > 100:
                vrt1 = "-".join(list(vt1))
                vrt2 = "-".join(list(vt2))
                fw.write(" " * 4 + "\"" + vrt1 + "\" -> \"" + vrt2 + "\"")
                fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=" + str(wgs[wg]) + "];\n")
            else:
                pass
                #fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=" + str(wgs[wg]) + " constraint=false];\n")


def print_dual_graph(db_cnt, cenid, out):
    mn_set = {x[0] for x in db_cnt.keys()}
    mn_set2 = {x[0] + "_" for x in db_cnt.keys()}

    cnt_mon = {x: 0 for x in mn_set}
    for x in db_cnt.keys():
        cnt_mon[x[0]] += db_cnt[x]
    cnt_mon2 = {x + "_": y for x, y in cnt_mon.items()}

    db_cnt2 = {(x[0], x[1] + "_"): y for x, y in db_cnt.items()}
    with open(os.path.join(out, cenid + ".dot"), "w") as fw:
        save_vert(fw, list(mn_set), list(mn_set2), cenid, cnt_mon, cnt_mon2)
        save_edges(fw, list(mn_set), list(mn_set2), db_cnt2)

    try:
        check_call(['dot', '-Tpng', os.path.join(out, cenid + ".dot"), '-o', os.path.join(out, cenid + ".png")])
    except Exception:
        return


def print_monomer_graph(db_cnt, cenid, matching, out):
    mn_set = {x[0] for x in db_cnt.keys()}

    cnt_mon = {x: 0 for x in mn_set}
    for x in db_cnt.keys():
        cnt_mon[x[0]] += db_cnt[x]


    with open(os.path.join(out, "mn_" + cenid + ".dot"), "w") as fw:
        save_vert_mn(fw, list(mn_set), cenid, cnt_mon)
        save_edges_mn(fw, list(mn_set), db_cnt, matching)
        fw.write("}\n")

    try:
        check_call(['dot', '-Tpng', os.path.join(out, "mn_" + cenid + ".dot"), '-o', os.path.join(out, "mn_" + cenid + ".png")])
    except Exception:
        return


def save_edges_sep_mn(fw, mns, sepdict):
    for vt1 in sorted(mns):
        if vt1 not in mns:
            continue
        vt1 = vt1[0]
        vser = vt1.split(".")[0]
        for v1, vt2, scr in sepdict:
            if v1 != vser:
                continue
            if (vt2,) not in mns:
                continue

            thr_wg = [1, 2, 5, 10, 500]
            wgs = [4, 3, 2, 1, 0]
            wg = 3
            while (scr < thr_wg[wg]):
                wg -= 1

            if scr > 10:
                continue

            fw.write(" " * 4 + "\"" + vt1 + "\" -> \"" + vt2 + "\"")
            fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=" + str(wgs[wg]) + " color=red constraint=true];\n")


def save_edges_posscore(fw, mns, posscore):
    for v1 in sorted(mns):
        for v2 in sorted(mns):
            vt1 = v1[0]
            vt2 = v2[0]
            if (vt1 == vt2):
                continue

            scr = posscore[(vt1, vt2)]

            if scr < 0.4:
                continue

            fw.write(" " * 4 + "\"" + vt1 + "\" -> \"" + vt2 + "\"")
            fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=2 color=orange constraint=true];\n")


def printk_graph(kcnt, cenid, matching, out, k, sepdict, posscore, thr=100):
    mn_set = {tuple(list(x)[:-1]) for x, y in kcnt.items() if y > thr} | \
             {tuple(list(x)[1:]) for x, y in kcnt.items() if y > thr}

    cnt_mon = {x: 0 for x in mn_set}
    for x in kcnt.keys():
        if tuple(list(x)[:-1]) in cnt_mon:
            cnt_mon[tuple(list(x)[:-1])] += kcnt[x]

    with open(os.path.join(out, "k" + str(k) + "_" + cenid + ".dot"), "w") as fw:
        save_vert_mn(fw, list(mn_set), cenid, cnt_mon)
        save_edges_mn(fw, list(mn_set), kcnt, matching)
        if k == 1:
            save_edges_sep_mn(fw, list(mn_set), sepdict)
            save_edges_posscore(fw, list(mn_set), posscore)
        fw.write("}\n")

    try:
        check_call(['dot', '-Tpng', os.path.join(out, "k" + str(k) + "_" + cenid + ".dot"), '-o',
                    os.path.join(out, "k" + str(k) + "_" + cenid + ".png")])
    except Exception:
        return

    pass


def print_k2_graph(trp_cnt, db_cnt, cenid, out):
    mn_set = {x[0] + "-" + x[1] for x,y in db_cnt.items() if y >= 100}
    cnt_mon = {x[0] + "-" + x[1] : y for x,y in db_cnt.items()}

    db2 = {(x[0] + "-" + x[1], x[1] + "-" + x[2]): y for x, y in trp_cnt.items()}
    with open(os.path.join(out, "k2_" + cenid + ".dot"), "w") as fw:
        save_vert_mn(fw, list(mn_set), cenid, cnt_mon)
        save_edges_mn(fw, list(mn_set), db2, {})

    try:
        check_call(['dot', '-Tpng', os.path.join(out, "k2_" + cenid + ".dot"), '-o',
                    os.path.join(out, "k2_" + cenid + ".png")])
    except Exception:
        return


def GetMaxMatching(kcnt):
    mn_set = {tuple(list(x)[:-1]) for x in kcnt.keys()} | {tuple(list(x)[1:]) for x in kcnt.keys()}

    cnt_mon = {x: 0 for x in mn_set}
    for x in kcnt.keys():
        cnt_mon[tuple(list(x)[:-1])] += kcnt[x]

    B = nx.Graph()
    B.add_nodes_from(list(mn_set), bipartite=0)
    B.add_nodes_from(list([(*x,"_") for x in mn_set]), bipartite=1)

    for x, w in kcnt.items():
        B.add_edge(tuple(list(x)[:-1]), tuple((list(x)[1:] + ["_"])), weight=-w)
    for x in mn_set:
        for y in mn_set:
            if list(x)[1:] != list(y)[:-1]:
                B.add_edge(x, (*y, "_"), weight=1000)
                continue

            if (*x, y[-1]) not in kcnt:
                B.add_edge(x, (*y, "_"), weight=0)

    matching = bipartite.matching.minimum_weight_full_matching(B, list(mn_set), "weight")
    return matching

def diffcenpos(cenv1, cenv2):
    mx1 = (0, 0)
    mx2 = (0, 0)

    for i in range(len(cenv1)):
        for j in range(len(cenv1[0])):
            if cenv1[i][j] > cenv1[mx1[0]][mx1[1]]:
                mx1 = (i, j)


    for i in range(len(cenv2)):
        for j in range(len(cenv2[0])):
            if cenv2[i][j] > cenv2[mx2[0]][mx2[1]]:
                mx2 = (i, j)

    return mx1[0] != mx2[0] and mx1[1] != mx2[1]

def mergeMon(sepdict, posscore, cenposvector):
    exch = {}
    for mn, mn2, cnt1 in sepdict:
        if (cnt1 < 6):
            if (mn, mn2) in posscore:
                #if diffcenpos(cenposvector[mn], cenposvector[mn2]):
                #    continue
                scr = posscore[(mn, mn2)]
                if scr > 0.3:
                    exch[mn] = mn2
    return exch


def getMnSim(mon):
    sim = {}
    for mn1 in mon:
        for mn2 in mon:
            sim[(mn1.id, mn2.id)] = min(SDutils.seq_identity(mn1.seq, mn2.seq),
                                        SDutils.seq_identity(mn1.seq, SDutils.rc(mn2.seq)))
    return sim


def handle_cen(cenid, args):
    maxk = 3
    k_cnt = calc_mn_order_stat(args.sdtsv, cenid, maxk=maxk)

    if args.sep != "-":
        df = pd.read_csv(os.path.join(args.sep, cenid + "all.csv")).values.tolist()
        sepdict = {("mn_" + str(x[1]), x[-1], x[-2]) for x in df}
    else:
        mons = SDutils.unique(SDutils.load_fasta(os.path.join(args.mon, cenid + "mn.fa")))
        sim = getMnSim(mons)
        sepdict = {(k[0], k[1], x) for k, x in sim.items() if x <= 10 and k[0] < k[1]}

    # PositionScore, cenvec = TriplesMatrix.handleAllMn(k_cnt[1], k_cnt[0])
    #
    # exch = mergeMon(sepdict, PositionScore, cenvec)
    #
    # k_cnt = calc_mn_order_stat(args.sdtsv, cenid, maxk=maxk, exchange=exch)
    # PositionScore, cenvec = TriplesMatrix.handleAllMn(k_cnt[1], k_cnt[0])
    # exch2 = mergeMon(sepdict, PositionScore, cenvec)
    # for x, y in exch2.items():
    #     exch[x] = y
    #
    # PrefixPosScore, cenvec = TriplesMatrix.PrefixPosScore(k_cnt[1], k_cnt[0])
    # exch3 = mergeMon(sepdict, PrefixPosScore, cenvec)
    # for x, y in exch3.items():
    #     if x not in exch:
    #         exch[x] = y
    #
    # SuffixPosScore, cenvec = TriplesMatrix.SuffixPosScore(k_cnt[1], k_cnt[0])
    # exch4 = mergeMon(sepdict, SuffixPosScore, cenvec)
    # for x, y in exch4.items():
    #     if x not in exch:
    #         exch[x] = y
    #
    # while True:
    #     update = False
    #     for x, y in exch.items():
    #         if y in exch:
    #             update = True
    #             exch[x] = exch[y]
    #     if not update:
    #         break
    #
    #
    # exchTrp = TriplesMatrix.SplitAllMn(k_cnt[1], k_cnt[0])
    # print("====EXCHANGE====")
    # print(exchTrp)

    #k_cnt = calc_mn_order_stat(args.sdtsv, cenid, maxk=maxk, exchange=exch, exchTrp=exchTrp)
    PositionScore, cenvec = TriplesMatrix.handleAllMn(k_cnt[1], k_cnt[0], thr=0)

    print(k_cnt)
    for k in range(1, maxk + 1):
        matching = GetMaxMatching(k_cnt[k - 1])
        #print_dual_graph(db_cnt, cenid, args.o)
        printk_graph(k_cnt[k - 1], cenid, matching, args.o, k, sepdict, PositionScore, thr=0)
        #print_monomer_graph(db_cnt, cenid, matching, args.o)
        #print_k2_graph(trp_cnt, db_cnt, cenid, args.o)


def main():
    args = parse_args()
    for i in range(4, 5):
        handle_cen("cen"+str(i)+"_", args)
    handle_cen("cenX_", args)


if __name__ == "__main__":
    main()