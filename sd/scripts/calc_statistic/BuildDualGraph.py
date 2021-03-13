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
    parser.add_argument("-o")

    return parser.parse_args()

def calc_mn_order_stat(sdtsv, cenid):
    db_cnt = {}
    with open(sdtsv, "r") as f:
        csv_reader = csv.reader(f, delimiter='\t')
        prev_row = []
        for row in csv_reader:
            if cenid not in row[0]:
                continue

            if row[2] == "start":
                continue

            if prev_row != []:
                pident = float(prev_row[4])
                pmon = prev_row[1]
                if pmon[-1] == "'":
                    pmon = pmon[:-1]
                identity = float(row[4])
                mon = row[1]
                if mon[-1] == "'":
                    mon = mon[:-1]

                if row[-1] != '?' and prev_row[-1] != '?':
                    if (pmon, mon) not in db_cnt:
                        db_cnt[(pmon, mon)] = 0
                    db_cnt[(pmon, mon)] += 1

            prev_row = row
    return db_cnt

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
            print(int(lg))
            fw.write(" " * 4 + '"' + vert + "\" [style=filled fillcolor=\"" + clr[int(lg)] + "\" label=\"" + vert + "[" + str(cnt_mn) + "]\"];\n")

    fw.write(" "*4 + "{rank = same; " + "; ".join(mn1) + ";}\n")
    fw.write(" "*4 + "{rank = same; " + "; ".join(mn2) + ";}\n")


def save_vert_mn(fw, mn1, cenName, cnt_mon1):
    fw.write("digraph " + cenName + " {\n")
    fw.write(" "* 4 + "rankdir=LR;\n")

    for mn, cnt_mon in [(mn1, cnt_mon1)]:
        for vert in mn:
            cnt_mn = 0
            lg = 0.01
            if vert in cnt_mon:
                cnt_mn = cnt_mon[vert]
                lg = math.log(cnt_mn)
            clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa"]
            print(int(lg))
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


def save_edges_mn(fw, mn1, db_cnt, matching):
    for vt1 in sorted(mn1):
        for vt2 in sorted(mn1):
            scr = 0
            if (vt1, vt2) in db_cnt:
                scr = db_cnt[(vt1, vt2)]

            thr_wg = [100000000, 1000, 500, 100, 1]
            wgs = [7, 5, 3, 1, 0]
            wg = 3
            while (scr > thr_wg[wg]):
                wg -= 1

            if matching[vt1] == vt2 + "_":
                fw.write(" " * 4 + "\"" + vt1 + "\" -> \"" + vt2 + "\"")
                fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=" + str(wgs[wg]) + " color=blue];\n")
            elif scr > 100:
                fw.write(" " * 4 + "\"" + vt1 + "\" -> \"" + vt2 + "\"")
                fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=" + str(wgs[wg]) + "];\n")
            else:
                pass
                #fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=" + str(wgs[wg]) + " constraint=false];\n")

    fw.write("}\n")


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

    try:
        check_call(['dot', '-Tpng', os.path.join(out, "mn_" + cenid + ".dot"), '-o', os.path.join(out, "mn_" + cenid + ".png")])
    except Exception:
        return


def GetMaxMatching(db_cnt):
    mn_set = {x[0] for x in db_cnt.keys()}

    cnt_mon = {x: 0 for x in mn_set}
    for x in db_cnt.keys():
        cnt_mon[x[0]] += db_cnt[x]

    B = nx.Graph()
    B.add_nodes_from(list(mn_set), bipartite=0)
    B.add_nodes_from(list([x + "_" for x in mn_set]), bipartite=1)

    for x, w in db_cnt.items():
        B.add_edge(x[0], x[1] + "_", weight=-w)
    for x in mn_set:
        for y in mn_set:
            if (x, y) not in db_cnt:
                B.add_edge(x, y + "_", weight=0)

    matching = bipartite.matching.minimum_weight_full_matching(B, list(mn_set), "weight")
    return matching



def handle_cen(cenid, args):
    db_cnt = calc_mn_order_stat(args.sdtsv, cenid)

    matching = GetMaxMatching(db_cnt)
    print_dual_graph(db_cnt, cenid, args.o)
    print_monomer_graph(db_cnt, cenid, matching, args.o)


def main():
    args = parse_args()
    for i in range(1, 23):
        handle_cen("cen"+str(i)+"_", args)
    handle_cen("cenX_", args)


if __name__ == "__main__":
    main()