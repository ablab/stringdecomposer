#!/usr/bin/env python3

import argparse
import edlib
import os
import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
import SDutils
import csv
import math
from matplotlib import pyplot as plt
from subprocess import check_call
import networkx as nx
from networkx.algorithms import bipartite

def parse_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-sdtsv")
    parser.add_argument("-sep")
    parser.add_argument("-o")

    return parser.parse_args()

def calc_mn_order_stat(sdtsv, cenid, maxk = 1):
    k_cnt = [{} for k in range(maxk)]
    RC_cnt = 0
    FR_cnt = 0
    rows = []
    with open(sdtsv, "r") as f:
        csv_reader = csv.reader(f, delimiter='\t')
        for row in csv_reader:
            if cenid not in row[0]:
                continue
            if row[2] == "start":
                continue

            rows.append(row)

    for i, row in enumerate(rows):
        cur_mons = ()
        identity = float(row[4])
        mon = row[1]
        if mon[-1] == "'":
            mon = mon[:-1]
            RC_cnt += 1
        else:
            FR_cnt += 1
        if row[-1] == '?':
            continue

        cur_mons = (mon,)

        for k in range(1, maxk + 1):
            if i - k >= 0 and rows[i - k] != []:
                pident = float(rows[i - k][4])
                pmon = rows[i - k][1]
                if pmon[-1] == "'":
                    pmon = pmon[:-1]

                if rows[i - k][-1] == '?':
                    break

                cur_mons = (*cur_mons, pmon)
                if cur_mons not in k_cnt[k - 1]:
                    k_cnt[k - 1][cur_mons] = 0
                k_cnt[k - 1][cur_mons] += 1
            else:
                break

    return k_cnt

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
            vertn = "-".join([x[3:] for x in list(vert)])
            vert = "-".join(list(vert))
            fw.write(" " * 4 + '"' + vert + "\" [style=filled fillcolor=\"" + clr[int(lg)] + "\" label=\"" + vertn + "[" + str(cnt_mn) + "]\"];\n")


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
        vt2 = sepdict[vt1][0]
        scr = sepdict[vt1][1]
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


def printk_graph(kcnt, cenid, matching, out, k, sepdict):
    mn_set = {tuple(list(x)[:-1]) for x, y in kcnt.items() if y > 100} | \
             {tuple(list(x)[1:]) for x, y in kcnt.items() if y > 100}

    cnt_mon = {x: 0 for x in mn_set}
    for x in kcnt.keys():
        if tuple(list(x)[:-1]) in cnt_mon:
            cnt_mon[tuple(list(x)[:-1])] += kcnt[x]

    with open(os.path.join(out, "k" + str(k) + "_" + cenid + ".dot"), "w") as fw:
        save_vert_mn(fw, list(mn_set), cenid, cnt_mon)
        save_edges_mn(fw, list(mn_set), kcnt, matching)
        if k == 1:
            save_edges_sep_mn(fw, list(mn_set), sepdict)
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


def handle_cen(cenid, args):
    maxk = 4
    k_cnt = calc_mn_order_stat(args.sdtsv, cenid, maxk=maxk)
    df = pd.read_csv(os.path.join(args.sep, cenid + ".csv")).values.tolist()
    sepdict = {"mn_" + str(x[1]) : (x[-1], x[-2]) for x in df}
    print(sepdict)

    for k in range(1, maxk + 1):
        matching = GetMaxMatching(k_cnt[k - 1])
        #print_dual_graph(db_cnt, cenid, args.o)
        printk_graph(k_cnt[k - 1], cenid, matching, args.o, k, sepdict)
        #print_monomer_graph(db_cnt, cenid, matching, args.o)
        #print_k2_graph(trp_cnt, db_cnt, cenid, args.o)


def main():
    args = parse_args()
    for i in range(1, 23):
        handle_cen("cen"+str(i)+"_", args)
    handle_cen("cenX_", args)


if __name__ == "__main__":
    main()