#!/usr/bin/env python3
import math

import networkx as nx
from networkx.drawing.nx_agraph import write_dot
from subprocess import check_call
import numpy as np
import SDutils

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


def canSplit(v, srunG, rG, epaths, k3cnt, handledV):
    oute = srunG.edges(v)
    ine = rG.edges(v)
    odeg = len(oute)
    ideg = len(ine)
    if ideg == 1 or odeg == 1:
        return True
    return False


def SplitV(v, mnrunG, rG, epaths, k3cnt, handledV, srunG):
    oute = list(srunG.edges(v))
    ine = list(rG.edges(v))
    odeg = len(oute)
    ideg = len(ine)

    alph = "".join([chr(ord('a') + i) for i in range(26)])
    vlist = []
    srunG.remove_node(v)
    rG.remove_node(v)
    for e0 in oute:
        for e1 in ine:
            vlist.append(v + alph[len(vlist)])
            srunG.add_node(vlist[-1])

            srunG.add_edge(e1[1], vlist[-1])
            print(v, e0[1])
            u = e0[1]
            if e0[1][-1] in alph:
                u = e0[1][:-1]

            srunG.add_edge(vlist[-1], e0[1])
            v1 = e1[1]
            if e1[1][-1] in alph:
                v1 = e1[1][:-1]

            u1 = v
            if len(ine) == 1:
                v1 = v
                u1 = u
            print(list(srunG.nodes()))
            print(list(srunG.edges()))
            print(srunG)
            print(e1)
            print(list(mnrunG.nodes()))
            srunG[e1[1]][vlist[-1]]["penwidth"] = mnrunG[v1][u1]["penwidth"]
            srunG[e1[1]][vlist[-1]]["label"] = mnrunG[v1][u1]["label"]
            srunG[vlist[-1]][e0[1]]["penwidth"] = mnrunG[v1][u1]["penwidth"]
            srunG[vlist[-1]][e0[1]]["label"] = mnrunG[v1][u1]["label"]

            rG.add_edge(vlist[-1], e1[1])
            rG.add_edge(e0[1], vlist[-1])
    return vlist



def SplitMnrunVert(mnrunG, epaths, k3cnt):
    srunG = nx.DiGraph()
    rG = nx.DiGraph()
    for v in mnrunG.nodes():
        rG.add_node(v)
        srunG.add_node(v)

    for e in mnrunG.edges():
        rG.add_edge(e[1], e[0])
        srunG.add_edge(e[0], e[1])
        srunG[e[0]][e[1]]["penwidth"] = mnrunG[e[0]][e[1]]["penwidth"]
        srunG[e[0]][e[1]]["label"] = mnrunG[e[0]][e[1]]["label"]

    handledV = {}

    for v in mnrunG.nodes():
        oute = mnrunG.edges(v)
        ine = rG.edges(v)
        odeg = len(srunG.edges(v))
        ideg = len(ine)

        if odeg < 2 and ideg < 2:
            handledV[v] = [v]
        elif canSplit(v, srunG, rG, epaths, k3cnt, handledV):
            handledV[v] = SplitV(v, mnrunG, rG, epaths, k3cnt, handledV, srunG)
        else:
            handledV[v] = [v]

    return srunG


def genCycleInner(mnrunG, prefixCycle, usedEdges, usedV, cycleList):
    def samecc(c1, c2):
        if len(c1) != len(c2):
            return False
        for j in range(len(c1) - 1):
            if c1[j:-1] + c1[:j] == c2[:-1]:
                return True
        return False

    if len(prefixCycle) > 2 and prefixCycle[-1] == prefixCycle[-2]:
        return

    if len(prefixCycle) > 1 and prefixCycle[0] == prefixCycle[-1]:
        for cc in cycleList:
            if samecc(cc, prefixCycle):
                break
        else:
            cycleList.append(prefixCycle)

    if len(prefixCycle) > 1 and prefixCycle[-1] == prefixCycle[-2]:
        return

    v = prefixCycle[-1]
    for e in mnrunG.edges(v):
        if e not in usedEdges and e[1] not in usedV:
            genCycleInner(mnrunG, prefixCycle + [e[1]], usedEdges | {e}, usedV, cycleList)


def genAllCycles(mnrunG):
    usedV = set()
    cycleList = []
    for v in mnrunG.nodes():
        if v not in usedV:
            genCycleInner(mnrunG, [v], set(), usedV, cycleList)
            usedV.add(v)

    return cycleList


def monomrunHOR2monomersHOR(cc, epaths):
    print(cc)
    mncc = []
    for i in range(len(cc) - 1):
        mncc += epaths[cc[i]][:-1]
    return mncc


def getHORcnt(mnHOR, monocen):
    cnt = 0
    for i in range(0, len(monocen) - len(mnHOR)):
        if mnHOR == monocen[i:i + len(mnHOR)]:
            cnt += 1
    return cnt


def isDecompose(cl, cls):
    isDecom = [False] * (len(cl))
    isDecom[len(cl) - 1] = True
    for j in range(len(cl) - 2, -1, -1):
        for k in range(len(cls)):
            if cls[k] == cl:
                continue
            if j + len(cls[k]) < len(cl) + 1:
                if cls[k] == cl[j:j + len(cls[k])]:
                    isDecom[j] = max(isDecom[j], isDecom[j + len(cls[k]) - 1])

    return isDecom[0]

def getHORs(exCls, cenid):
    hors = []
    for i in range(len(exCls)):
        if not isDecompose(exCls[i], exCls):
            nm = "C" + cenid[3:-1] + "H" + str(len(hors) + 1)
            hors.append((nm, exCls[i]))

    return hors


def getHORcntMR(HOR, mnrunG):
    mnEdge = 100000
    for i in range(len(HOR) - 1):
        mnEdge = min(mnEdge, int(mnrunG[HOR[i]][HOR[i + 1]]["label"]))
    return mnEdge


def CAtoIA(CAIA, epath):
    return [CAIA[mn] for mn in epath]


def BuildAndShowMonorunGraph(k2cnt, k3cnt, ofile, monocen, cenid, CAIA, vLim=100, eLim = 100):
    vcnt = {v : 0 for v, u in k2cnt.keys()}
    for vu, cnt in k2cnt.items():
        vcnt[vu[0]] += cnt

    #print(vcnt)

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

    def smpl(vr):
        plv = vr.split("+")
        rplv = []
        for cv in plv:
            rrv = ""
            while len(rrv) + 1 < len(cv) and cv[:len(rrv) + 1].isalpha():
                rrv += cv[len(rrv)]

            if len(cv.split(".")) > 1:
                rrv += cv.split(".")[-1]
            rplv.append(rrv)

        rv = "+".join(rplv)
        if len(rv) > 1:
            rv = f'({rv})'
        return rv

    epaths = {}
    for le in lesall:
        with open("L.csv", "a") as fw:
            tmpia = CAtoIA(CAIA, le.epath)
            fw.write(le.name + "\t" + "".join([smpl(vr) for vr in le.epath]) + "\t" + str(tmpia[0]) + "," + ",".join([vl.split('.')[-1] for vl in tmpia[1:]]) + "\n")
        epaths[le.name] = le.epath
        for le2 in lesall:
            if le.epath[-1] == le2.epath[0]:
                #print(len(le.epath), len(le2.epath))

                #print(le2.epath)
                if (le.epath[-2], le2.epath[0], le2.epath[1]) in k3cnt and \
                        k3cnt[(le.epath[-2], le2.epath[0], le2.epath[1])] >= eLim:
                    mnrunG.add_edge(le.name, le2.name,
                                    penwidth=1.5*(math.log(float(k3cnt[(le.epath[-2], le2.epath[0], le2.epath[1])])) - 4),
                                    label=str(k3cnt[(le.epath[-2], le2.epath[0], le2.epath[1])]))

    write_dot(mnrunG, ofile)
    try:
        check_call(['dot', '-Tpng', ofile, '-o', ".".join(ofile.split(".")[:-1]) + ".png"])
    except Exception:
        return

    srunG = SplitMnrunVert(mnrunG, epaths, k3cnt)
    sruno = ".".join(ofile.split(".")[:-1]) + "_splv.dot"
    srunopng = ".".join(ofile.split(".")[:-1]) + "_splv.png"
    write_dot(srunG, sruno)
    try:
        check_call(['dot', '-Tpng', sruno, '-o', srunopng])
    except Exception:
        return



    clsList = genAllCycles(mnrunG)
    mnPotentialHOR = [monomrunHOR2monomersHOR(cc, epaths) for cc in clsList]
    HORcnt = [getHORcnt(mnHOR + mnHOR, monocen) for mnHOR in mnPotentialHOR]
    print("AllCycles:")
    print(np.array(clsList))
    print("potential HOR:")
    with open("HORs.tsv", "a") as fw:
        for i in range(len(HORcnt)):
            fw.write(str(clsList[i]) +  "\t" + str(mnPotentialHOR[i]) +  "\t" + str(HORcnt[i]) + "\n")

    exCls = [clsList[i] for i in range(len(clsList)) if HORcnt[i] > 0]
    HORs = getHORs(exCls, cenid)
    print(HORs)
    with open("HORs/HOR" + cenid[3:-1] + ".tsv", "w") as fw:
        for i in range(len(HORs)):
            hor = monomrunHOR2monomersHOR(HORs[i][1], epaths)
            fw.write(HORs[i][0] + "\t" + ",".join(hor) + "\t" + ",".join(HORs[i][1]) + "\n")

    with open("HORscnt.tsv", "a") as fw:
        for i in range(len(HORs)):
            hor = monomrunHOR2monomersHOR(HORs[i][1], epaths)
            fw.write(HORs[i][0] + ": "+"+".join(HORs[i][1][:-1]) + "(" + str(len(hor)) + "-mer)" + "\t" + str(getHORcntMR(HORs[i][1], mnrunG)) + "\n")