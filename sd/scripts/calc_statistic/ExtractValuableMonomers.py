#!/usr/bin/env python3

import argparse
import edlib
import os
import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
from SDutils import get_blocks
from SDutils import seq_identity
from SDutils import  rc

def get_sq_sum(blocks, dists):
    dists = np.array(dists)
    dists = dists.min(axis=1)
    dists **=2
    return (sum(dists)/len(blocks))**0.5


def get_subset_sq_sum(blocks, dists, mns, used_mn):
    if len(used_mn) == 0:
        return 100
    sub_dists = [[dists[j][i] for i, mn in enumerate(mns) if mn.id in used_mn] for j in range(0, len(blocks))]
    return get_sq_sum(blocks, sub_dists)


def getValuableMonomers(path_seq, tsv_res, CAmn, odir):
    blocks = get_blocks(path_seq, tsv_res)
    dists = [[min(seq_identity(str(mn.seq), bl), seq_identity(str(mn.seq), rc(bl))) for mn in CAmn] for bl in blocks]

    mnCnt = [0] * len(CAmn)
    for i in range(len(blocks)):
        mnCnt[dists[i].index(min(dists[i]))] += 1

    OrderMn = list(zip(mnCnt, [i for i in range(len(CAmn))]))
    OrderMn.sort(reverse=True)

    mnStatistic = []

    used_mn = set()
    all_mn = set()
    cfreq = set()
    CntAll = 0
    CntVal = 0

    pr_sqsum = 100
    lstFrq = 0
    for i in range(0, len(CAmn)):
        curMn = CAmn[OrderMn[i][1]].id
        idmn = OrderMn[i][1]

        mnStatistic.append([curMn, mnCnt[idmn]/len(blocks), (CntAll + mnCnt[idmn])/len(blocks), CntVal/len(blocks), "distortAll", "distorVal", False])
        sq_sumUsed = get_subset_sq_sum(blocks, dists, CAmn, used_mn | {curMn})
        sq_allUsed = get_subset_sq_sum(blocks, dists, CAmn, all_mn | {curMn})
        mnStatistic[-1][-3] = sq_allUsed/get_subset_sq_sum(blocks, dists, CAmn, all_mn)
        mnStatistic[-1][-2] = sq_sumUsed/get_subset_sq_sum(blocks, dists, CAmn, used_mn)
        if CntAll/len(blocks) < 0.90 or lstFrq < 1.5*mnCnt[idmn]/len(blocks):
            cfreq.add(curMn)
            lstFrq = mnCnt[idmn]/len(blocks)

        CntAll += mnCnt[idmn]
        all_mn.add(curMn)


        if pr_sqsum*0.99 > sq_sumUsed:
            used_mn.add(curMn)
            pr_sqsum = sq_sumUsed
            mnStatistic[-1][-1] = True
            mnStatistic[-1][-4] = (CntVal + mnCnt[idmn]) / len(blocks)
            CntVal += mnCnt[idmn]

    df = pd.DataFrame(mnStatistic, columns=["MonomerID", "Frequence", "Cumulative frequences",
                                            "Valuable Cumulative frequences", "ALLdistortionRatio",
                                            "ValueDistortionRatio", "IsValuable"])
    df.to_csv(os.path.join(odir, "MonomerStat.csv"))

    MonomersNew = [mn for mn in CAmn if mn.id in cfreq]
    return MonomersNew