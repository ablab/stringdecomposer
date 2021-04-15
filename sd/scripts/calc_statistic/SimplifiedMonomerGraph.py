#!/bin/usr/env python3
import networkx as nx
from networkx.algorithms import bipartite
from networkx.drawing.nx_agraph import write_dot
from subprocess import check_call
import math
import os


def GetMaxMatching(kcnt):
    mn_set = {tuple(list(x)[:-1]) for x in kcnt.keys()} | {tuple(list(x)[1:]) for x in kcnt.keys()}

    cnt_mon = {x: 0 for x in mn_set}
    for x in kcnt.keys():
        cnt_mon[tuple(list(x)[:-1])] += kcnt[x]

    B = nx.Graph()
    B.add_nodes_from(list(mn_set), bipartite=0)
    B.add_nodes_from(list([(*x,"_") for x in mn_set]), bipartite=1)

    print("nodes", B.nodes())
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


def addNode(G, mn, mncnt, IAnm, IsHybrid):
    lg = 0.01
    if mncnt > 0:
        lg = math.log(mncnt)
    clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa"]
    curc = "pink" if IsHybrid else clr[int(lg)]
    G.add_node(mn, style=f'filled', fillcolor=f'{curc}', label=f'{mn}({IAnm}) [{str(int(mncnt))}]')

def addEdges(G, mn1, mn2, cnt, thr):
    scr = cnt
    thr_wg = [100000000, 1000, 500, 100, 1]
    wgs = [7, 5, 3, 1, 0]
    wg = 3
    while (scr > thr_wg[wg]):
        wg -= 1

    if scr > thr:
        G.add_edge(mn1, mn2, label=f'{scr}', penwidth=f'{wgs[wg]}')


def PrintSimplifiedGraph(kcnt, mncnt, CAIA, HybridSet, outd, cenid, edgeThr=100):
    matching = GetMaxMatching(kcnt)
    G = nx.DiGraph()
    for mn in mncnt.keys():
        addNode(G, mn, mncnt[mn], CAIA[mn], False)

    print(matching)
    cycleid = {mn: i for i, mn in enumerate(list(mncnt.keys()))}
    def redrw(mn1, mn2):
        c1 = cycleid[mn1]
        c2 = cycleid[mn2]

        mnclr = min(c1, c2)
        for mn in cycleid:
            if cycleid[mn] == c1 or cycleid[mn] == c2:
               cycleid[mn] = mnclr


    for mn1, mn2 in matching.items():
        if len(mn1) > 1:
            continue

        if (mn1[0], mn2[0]) in kcnt:
            redrw(mn1[0], mn2[0])
            addEdges(G, mn1[0], mn2[0], kcnt[(mn1[0], mn2[0])], edgeThr)

    print("CycleID", cycleid)
    for mn1, mn2 in kcnt.keys():
        if cycleid[mn1] != cycleid[mn2]:
            addEdges(G, mn1, mn2, kcnt[mn1, mn2], edgeThr)

    sruno = os.path.join(outd, "simpl" + cenid + ".dot")
    srunopng = os.path.join(outd, "simpl" + cenid + ".png")
    write_dot(G, sruno)
    try:
        check_call(['dot', '-Tpng', sruno, '-o', srunopng])
    except Exception:
        return
    pass