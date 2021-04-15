#!/usr/bin/env python3
import SDutils

def get_hybrid_len(main_mn, mn1, mn2):
    resDiv = 500
    mn_identity = 100
    bst_res = (0, 0)
    for prfx in range(30, len(mn1.seq)):
        suffix = len(str(main_mn.seq)) - prfx
        if suffix < 30:
            break
        hbr = str(mn1.seq)[:prfx] + str(mn2.seq)[len(mn2.seq) - suffix:]
        cur_identity = SDutils.seq_identity(hbr, str(main_mn.seq))
        if mn_identity > cur_identity:
            mn_identity = cur_identity
            bst_res =  (prfx, suffix)

    if (mn_identity * 2 <= resDiv):
        return (bst_res[0], bst_res[1], mn_identity)
    return (0, 0, 100)


def isHybrid(main_mn, mn1, mn2):
    pr, sf, idn = get_hybrid_len(main_mn, mn1, mn2)
    if idn < 6:
        return True
    return False

def getHybridINFO(mnpath, vcnt):
    hybridSet = set()
    monCA = SDutils.load_fasta(mnpath)
    for i in range(len(monCA)):
        for j in range(len(monCA)):
            for g in range(len(monCA)):
                if i != j and j != g and i != g:
                    if vcnt[monCA[j].id] > vcnt[monCA[i].id] and vcnt[monCA[g].id] > vcnt[monCA[i].id]:
                        if isHybrid(monCA[i], monCA[j], monCA[g]):
                            hybridSet.add(monCA[i].id)
    return hybridSet
