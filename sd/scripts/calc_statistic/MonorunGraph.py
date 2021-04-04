#!/usr/bin/env python3
import math

import networkx as nx
from networkx.drawing.nx_agraph import write_dot
from subprocess import check_call

class LongEdge:
    def __init__(self):
        self.epath = []
        self.weight = 0
        self.name = ""

def initLongEdgeWeight(outLongE, k2cnt):
    return


def initLongEdgeNames(outLongE):
    namecnt = {}
    for v, les in outLongE.items():
        for le in les:
            cname = "L" + str(len(le.epath) - 1)
            if le.epath[0] == le.epath[-1]:
                cname = "C" + str(len(le.epath) - 1)
            if cname not in namecnt:
                namecnt[cname] = 1
                le.name = cname
            else:
                le.name = cname + "-" + str(namecnt[cname])
                namecnt[cname] += 1


def BuildAndShowMonorunGraph(k2cnt, k3cnt, ofile, vLim=100, eLim = 100):
    vcnt = {v : 0 for v, u in k2cnt.keys()}
    for vu, cnt in k2cnt.items():
        vcnt[vu[0]] += cnt

    print(vcnt)

    ine = {v : [] for v, cnt in vcnt.items() if cnt >= vLim}
    oute = {v : [] for v, cnt in vcnt.items() if cnt >= vLim}

    for vu, cnt in k2cnt.items():
        if cnt < eLim:
            continue
        ine[vu[1]].append(vu[0])
        oute[vu[0]].append(vu[1])

    usedv = set()

    outLongE = {v: [] for v in oute.keys()}
    for v in ine.keys():
        if v not in usedv and (len(oute[v]) != 1 or len(ine[v]) != 1):
            usedv |= {v}
            for u in oute[v]:
                curE = LongEdge()
                outLongE[v].append(curE)
                curE.epath.append(v)
                cu = u
                while len(oute[cu]) == 1 and len(ine[cu]) == 1:
                    usedv.add(cu)
                    curE.epath.append(cu)
                    cu = oute[cu][0]
                curE.epath.append(cu)

    for v in ine.keys():
        if v not in usedv:
            usedv |= {v}
            curE = LongEdge()
            outLongE[v].append(curE)
            curE.epath.append(v)
            u = oute[v][0]
            while u != v:
                usedv |= {u}
                curE.epath.append(u)
                u = oute[u][0]
            curE.epath.append(v)

    initLongEdgeWeight(outLongE, k2cnt)
    initLongEdgeNames(outLongE)


    lesall = []
    mnrunG = nx.DiGraph()
    for v, les in outLongE.items():
        for le in les:
            lesall.append(le)
            mnrunG.add_node(le.name)

    for le in lesall:
        print(le.epath)
        for le2 in lesall:
            if le.epath[-1] == le2.epath[0]:
                print(len(le.epath), len(le2.epath))

                print(le2.epath)
                if k3cnt[(le.epath[-2], le2.epath[0], le2.epath[1])] >= eLim:
                    mnrunG.add_edge(le.name, le2.name,
                                    penwidth=1.5*(math.log(float(k3cnt[(le.epath[-2], le2.epath[0], le2.epath[1])])) - 4),
                                    label=str(k3cnt[(le.epath[-2], le2.epath[0], le2.epath[1])]))

    write_dot(mnrunG, ofile)
    try:
        check_call(['dot', '-Tpng', ofile, '-o', ".".join(ofile.split("."))[:-1] + ".png"])
    except Exception:
        return
