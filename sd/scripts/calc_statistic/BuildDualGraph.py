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

from MonorunGraph import BuildAndShowMonorunGraph

import TriplesMatrix
from TriplesMatrix import calc_mn_order_stat
import SDutils
import HybridUtils
import SimplifiedMonomerGraph

def parse_args():
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-sdtsv")
    parser.add_argument("-seq")
    parser.add_argument("-sep", default="-")
    parser.add_argument("-mon", default="-")
    parser.add_argument("-monIA", default="-")
    parser.add_argument("--blue", dest="blue", action='store_true')
    parser.add_argument("--red", dest="red", action="store_true")
    parser.add_argument("--norm", dest="norm", action="store_true")
    parser.add_argument("--maxk", dest="maxk", default=1, type=int)
    parser.add_argument("--monorun", dest="monorun", action="store_true")
    parser.add_argument("-o")

    return parser.parse_args()


def save_vert_mn(fw, mn1, cenName, cnt_mon1, CAIA, HybridSet):
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
            clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa", "#f5fffa", "#f5fffa", "#f5fffa", "#f5fffa"]
            print(int(lg))
            #vertn = "-".join([x[3:] for x in list(vert)])
            vert = "-".join(list(vert))
            curc = clr[int(lg)]
            #if vert in HybridSet:
            #    curc = "pink"
            if vert in CAIA:
                fw.write(f'    "{vert}" [style=filled fillcolor="{curc}" label="{vert}({CAIA[vert]}) [{str(int(cnt_mn))}]"];\n')
            else:
                fw.write(f'    "{vert}" [style=filled fillcolor="{curc}" label="{vert}[{str(int(cnt_mn))}]"];\n')


def save_edges_mn(fw, mn1, kcnt, matching, thr=100):
    ecnt = 0
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

            if scr > thr:
                ecnt += 1

            if matching.get(vt1, "") == (*vt2, "_"):
                vrt1 = "-".join(list(vt1))
                vrt2 = "-".join(list(vt2))
                fw.write(" " * 4 + "\"" + vrt1 + "\" -> \"" + vrt2 + "\"")
                fw.write(" [label=\"" +  str(int(scr)) + "\" penwidth=" + str(wgs[wg]) + " color=blue];\n")
            elif scr > thr:
                vrt1 = "-".join(list(vt1))
                vrt2 = "-".join(list(vt2))
                fw.write(" " * 4 + "\"" + vrt1 + "\" -> \"" + vrt2 + "\"")
                fw.write(" [label=\"" + str(int(scr)) + "\" penwidth=" + str(wgs[wg]))
                if scr < 5:
                    fw.write(" constraint = false];\n")
                else:
                    fw.write("];\n")
            else:
                pass
                #fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth=" + str(wgs[wg]) + " constraint=false];\n")
    return ecnt


def save_edges_sep_mn(fw, mns, sepdict):
    print("EdgeSep")
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

            thr_wg = [-1, 1, 2, 5, 11, 5000]
            wgs = [5, 4, 3, 2, 1, 0]
            wg = 4

            if scr > 10:
                continue

            while (scr < thr_wg[wg]):
                wg -= 1

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


def printk_graph(kcnt, cnt_mon, cenid, matching, out, k, sepdict, posscore, args, CAIA, HybridINFO, thr=100, edgeThr=100):
    mn_set = {tuple(list(x)[:-1]) for x, y in kcnt.items() if y > thr} | \
             {tuple(list(x)[1:]) for x, y in kcnt.items() if y > thr}

    if k == 1:
        cntmn2 = {(mn,): cnt for mn, cnt in cnt_mon.items()}
        cnt_mon = cntmn2

    vcnt = len(mn_set)
    ecnt = 0
    with open(os.path.join(out, "k" + str(k) + "_" + cenid + ".dot"), "w") as fw:
        save_vert_mn(fw, list(mn_set), cenid, cnt_mon, CAIA, HybridINFO)
        ecnt += save_edges_mn(fw, list(mn_set), kcnt, matching, edgeThr)
        if k == 1 and args.red:
             save_edges_sep_mn(fw, list(mn_set), sepdict)
             save_edges_posscore(fw, list(mn_set), posscore)
        fw.write("}\n")

    with open(os.path.join(out, "MG_vecnt.csv"), "a") as fw:
        fw.write(f'{cenid}, {vcnt}, {ecnt}\n')

    try:
        check_call(['dot', '-Tpng', os.path.join(out, "k" + str(k) + "_" + cenid + ".dot"), '-o',
                    os.path.join(out, "k" + str(k) + "_" + cenid + ".png")])
    except Exception:
        return


def getMnSim(mon):
    sim = {}
    for mn1 in mon:
        for mn2 in mon:
            sim[(mn1.id, mn2.id)] = min(SDutils.seq_identity(mn1.seq, mn2.seq),
                                        SDutils.seq_identity(mn1.seq, SDutils.rc(mn2.seq)))
    return sim


def getEdheThr(k2cnt):
    vcnt = {v : 0 for v, u in k2cnt.keys()}
    for vu, cnt in k2cnt.items():
        vcnt[vu[0]] += cnt

    print(vcnt)

    minW = min(100, min([cnt for v, cnt in vcnt.items()])*0.9)
    return minW


def mapCAIAmn(mnCA, mnIA):
    monCA = SDutils.load_fasta(mnCA)
    monIA = SDutils.load_fasta(mnIA)

    mapCAIA = {}
    for mCA in monCA:
        bstidn = 100
        for mIA in monIA:
            bstidn = min(SDutils.seq_identity(str(mCA.seq), str(mIA.seq)), SDutils.seq_identity(str(mCA.seq), SDutils.rc(mIA.seq)), bstidn)

        for mIA in monIA:
            if SDutils.seq_identity(str(mCA.seq), str(mIA.seq)) == bstidn:
                mapCAIA[mCA.id] = mIA.id
            if SDutils.seq_identity(str(mCA.seq), SDutils.rc(mIA.seq)) == bstidn:
                mapCAIA[mCA.id] = mIA.id + "-rev"
    return mapCAIA

def normalize(k_cnt):
    maxCnt = max(k_cnt[0].values())
    for i in range(len(k_cnt)):
        for ky in k_cnt[i].keys():
             k_cnt[i][ky] = k_cnt[i][ky]/maxCnt*100
    return k_cnt


def handle_cen(cenid, args):
    maxk = args.maxk
    k_cnt = calc_mn_order_stat(os.path.join(args.sdtsv, cenid + "dec.tsv"), cenid, maxk=max(2, maxk))
    if args.norm:
        k_cnt = normalize(k_cnt)

    edgeThr = getEdheThr(k_cnt[1])

    if args.sep != "-":
        df = pd.read_csv(os.path.join(args.sep, cenid + "all.csv")).values.tolist()
        sepdict = {("mn_" + str(x[1]), x[-1], x[-2]) for x in df}
    else:
        mons = SDutils.unique(SDutils.load_fasta(os.path.join(args.mon, cenid + "mn.fa")))
        sim = getMnSim(mons)
        sepdict = {(k[0], k[1], x) for k, x in sim.items() if x <= 10 and k[0] < k[1]}

    CAIA = mapCAIAmn(os.path.join(args.mon, cenid + "mn.fa"), os.path.join(args.monIA, cenid + "mn.fa"))
    with open("map_" + cenid[:-1] + ".tsv", "w") as fw:
        for ca, ia in CAIA.items():
            fw.write(f'{ca}\t{ia}\n')
    HybridINFO = HybridUtils.getHybridINFO(os.path.join(args.mon, cenid + "mn.fa"), k_cnt[0])
    print("HybridSet:", HybridINFO)

    PositionScore, cenvec = TriplesMatrix.handleAllMn(k_cnt[2], k_cnt[1], thr=0)

    for k in range(1, maxk + 1):
        matching = {}
        if args.blue:
            matching = SimplifiedMonomerGraph.GetMaxMatching(k_cnt[k])
        printk_graph(k_cnt[k], k_cnt[k - 1], cenid, matching, args.o, k, sepdict, PositionScore, args, CAIA, HybridINFO, thr=0, edgeThr=edgeThr)

    mncen = SDutils.get_monocent(os.path.join(args.sdtsv, cenid + "dec.tsv"))
    if args.monorun:
        BuildAndShowMonorunGraph(k_cnt[1], k_cnt[2], os.path.join(args.o, cenid + "mnrun.dot"), mncen, cenid, CAIA, vLim=0, eLim=edgeThr)
        SimplifiedMonomerGraph.PrintSimplifiedGraph(k_cnt[1], k_cnt[0], CAIA, HybridINFO, args.o, cenid,
                                                    os.path.join(args.sdtsv, cenid + "dec.tsv"),
                                                    os.path.join(args.seq, cenid[:-1] + "ct.fa"),
                                                    os.path.join(args.mon, cenid + "mn.fa"),
                                                    edgeThr=edgeThr)

def main():
    args = parse_args()
    for i in range(1, 23):
        try:
            handle_cen("cen"+str(i)+"_", args)
        except Exception:
            continue

    handle_cen("cenX_", args)


if __name__ == "__main__":
    main()