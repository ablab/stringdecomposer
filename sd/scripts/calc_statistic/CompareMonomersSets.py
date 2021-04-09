#!/usr/bin/env python3

import argparse
import edlib
import os
import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
import SDutils
from SDutils import get_blocks
from SDutils import  rc
from SDutils import sys_call
from SDutils import unique
from SDutils import load_fasta
from SDutils import run_SD
from SDutils import seq_identity
import ExtractValuableMonomers

res_str = ""

pathToSD = "/Bmo/kolga/soft/stringdecomposer/sd/run_decomposer.py"


def local_dist(a, b):
    result = edlib.align(b, a, mode="HW", task="locations")
    if result["editDistance"] == -1:
        return 10**9
    return result["editDistance"]*100/len(b)


def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--CA", dest="camon")
    parser.add_argument("--IA", dest="ivanmon")
    parser.add_argument("--seq", dest="seq")
    parser.add_argument("-o", dest = "outdir")

    return parser.parse_args()


def getSameMons(CA_mon, Ivan_mon):
    Bsubset = []
    for mn in Ivan_mon:
        dst_list = [min(seq_identity(str(x.seq), str(mn.seq)), seq_identity(rc(str(x.seq)), str(mn.seq))) for x in CA_mon]
        if min(dst_list) < 5:
            Bsubset.append(dst_list.index(min(dst_list)))

    Bsubset = list(set(Bsubset))
    return Bsubset


def run_SD_mn(mnList, seqPath, outName):
    os.makedirs(outName, exist_ok=True)
    origin_dir = os.path.abspath(os.getcwd())
    os.chdir(outName)
    if not os.path.exists("final_decomposition.tsv"):
        pathToMon = "mon.fa"
        with open(pathToMon, "w") as fw:
            for x in mnList:
                SeqIO.write(x, fw, "fasta")
        sys_call(["python3", pathToSD, seqPath, pathToMon, "-t", "30", "--fast"])
    os.chdir(origin_dir)
    return os.path.join(outName, "final_decomposition.tsv")


def addMostFreq(CA_mon, BsubsetIDX, tsv_res, cntAdd):
    df_sd = pd.read_csv(tsv_res, "\t")
    mn_count = {(len(df_sd.loc[df_sd.iloc[:, 1] == mn.id]), mn.id) for mn in CA_mon}
    BsubsetName = {CA_mon[x].id for x in BsubsetIDX}
    cntAdd += len(BsubsetName)
    for cnt, mn_nm in sorted(mn_count)[::-1]:
        if len(BsubsetName) < cntAdd:
            BsubsetName.add(mn_nm)
    resMN = [x for x in CA_mon if x.id in BsubsetName]
    return resMN


def get_squere_mean(tsv_path, mns):
    df_sd = pd.read_csv(tsv_path, "\t")
    return (sum((np.array([100] * len(df_sd)) - np.array(df_sd.iloc[:, 4]))**2)/len(df_sd))**(1/2)


def get_blocks_interseq(path_seq, tsva, tsvb):
    seqs_dict = {}
    for record in SeqIO.parse(path_seq, "fasta"):
        seqs_dict[record.id] = str(record.seq).upper()

    df1 = pd.read_csv(tsva, "\t", header=None)
    df2 = pd.read_csv(tsvb, "\t", header=None)

    blocks_inter = []
    i = 0
    j = 0
    while i < len(df1) and j < len(df2):
        if float(df1.iloc[i, 4]) < 60:
            i += 1
            continue

        if float(df2.iloc[j, 4]) < 60:
            j += 1
            continue

        lft = max(df1.iloc[i, 2], df2.iloc[j, 2])
        rgt = min(df1.iloc[i, 3], df2.iloc[j, 3])

        if rgt + 1 - lft > 150:
            blocks_inter.append(seqs_dict[df1.iloc[i, 0]][lft:rgt + 1])

        if df1.iloc[i, 2] < df2.iloc[j, 2]:
            i += 1
        else:
            j += 1

    return blocks_inter


def get_sq_sum(blocks, dists):
    dists = np.array(dists)
    dists = dists.min(axis=1)
    dists **=2
    return (sum(dists)/len(blocks))**0.5


def get_subset_sq_sum(blocks, dists, mns, used_mn):
    sub_dists = [[dists[j][i] for i, mn in enumerate(mns) if mn.id in used_mn] for j in range(0, len(blocks))]
    return get_sq_sum(blocks, sub_dists)


def update_used_mn(blocks, dists, CAmn, used_mn):
    bst_mn = -1
    bst_sc = 100
    for mn in CAmn:
        if mn.id not in used_mn:
            sd = get_subset_sq_sum(blocks, dists, CAmn, used_mn | set({mn.id}))
            if sd < bst_sc:
                bst_sc= sd
                bst_mn = mn.id

    used_mn.add(bst_mn)
    return bst_mn


def elbow_calc(path_seq, tsv_B_res, Bmn, CAmn, f):
    blocks = get_blocks(path_seq, tsv_B_res)

    dists = [[min(seq_identity(str(mn.seq), bl), seq_identity(str(mn.seq), rc(bl))) for mn in CAmn] for bl in blocks]

    used_mn = { mn.id for mn in Bmn }

    elbow_str = ""
    prev_mn = ""
    for i in range(len(Bmn), len(CAmn) + 1):
        if i > 0:
            sq_sum = get_subset_sq_sum(blocks, dists, CAmn, used_mn)
            print("{}({}):\t{:.4f};\n".format(i, prev_mn, sq_sum))
            elbow_str += "{}({}):\t{:.4f};\n".format(i, prev_mn, sq_sum)

        if i < len(CAmn):
            prev_mn = update_used_mn(blocks, dists, CAmn, used_mn)

    f.write("elbow " + elbow_str + "\n")


def getValuableMonomers(path_seq, tsv_res, CAmn):
    blocks = get_blocks(path_seq, tsv_res)
    dists = [[min(seq_identity(str(mn.seq), bl), seq_identity(str(mn.seq), rc(bl))) for mn in CAmn] for bl in blocks]

    used_mn = set()

    prev_mn = ""
    pr_sqsum = 100
    for i in range(0, len(CAmn) + 1):
        if i > 0:
            sq_sum = get_subset_sq_sum(blocks, dists, CAmn, used_mn)
            if pr_sqsum < sq_sum * 1.05:
                used_mn.remove(prev_mn)
                MonomersNew = [mn for mn in CAmn if mn.id in used_mn]
                return MonomersNew
            pr_sqsum = sq_sum

        if i < len(CAmn):
            prev_mn = update_used_mn(blocks, dists, CAmn, used_mn)

    MonomersNew = [mn for mn in CAmn if mn.id in used_mn]
    print(MonomersNew)
    return MonomersNew


def calc_mean_for_subseq(seq, tsv_B_res, tsv_Ivan_res, CAmn, IAmn, f):
    blocks = get_blocks_interseq(seq, tsv_B_res, tsv_Ivan_res)
    f.write("#blocks after interseq: " + str(len(blocks)) + "\n")

    dists_CA = [[min(local_dist(str(mn.seq), bl), local_dist(str(mn.seq), rc(bl))) for mn in CAmn] for bl in blocks]
    dists_IA = [[min(local_dist(str(mn.seq), bl), local_dist(str(mn.seq), rc(bl))) for mn in IAmn] for bl in blocks]

    CA_sq_sum = get_sq_sum(blocks, dists_CA)
    IA_sq_sum = get_sq_sum(blocks, dists_IA)

    f.write("CA sq sum for interseq blocks: " + str(CA_sq_sum) + "\n")
    f.write("IA sq sum for interseq blocks: " + str(IA_sq_sum) + "\n")

    global res_str
    res_str += "\t" + "{:.4f}".format(IA_sq_sum)
    res_str += "\t" + "{:.4f}".format(CA_sq_sum)

    return CA_sq_sum, IA_sq_sum

def blocks_stats(res_tsv1, res_tsv2, f):
    df1 = pd.read_csv(res_tsv1, "\t", header=None)
    df2 = pd.read_csv(res_tsv2, "\t", header=None)

    f.write("CA blocks#: " + str(len(df1)) + "\n")
    f.write("IA blocks#: " + str(len(df2)) + "\n")

    df3 = pd.merge(df1, df2, how="inner", on=[0,2,3])
    f.write("SAME blocks#: " + str(len(df3)) + "\n")

    CAmns = set(df1.iloc[:, 1])
    IAmns = set(df2.iloc[:, 1])

    for mn in CAmns:
        subdf = df1.loc[df1.iloc[:, 1] == mn]
        f.write("CA blocks# for " + mn + ": " + str(len(subdf)) + "\n")
        f.write("SAME CA and IA blocks # for " + mn + ": " + str(len(pd.merge(subdf, df2, how="inner", on=[0,2,3]))) + "\n")

    for mn in IAmns:
        subdf = df2.loc[df2.iloc[:, 1] == mn]
        f.write("IA blocks# for " + mn + ": " + str(len(subdf)) + "\n")
        f.write("SAME CA and IA blocks # for " + mn + ": " + str(len(pd.merge(subdf, df1, how="inner", on=[0,2,3]))) + "\n")

    f.write("\n\n")


def blocks_sqmean(path_seq, tsv, mns, f, name):
    blocks = get_blocks(path_seq, tsv)
    f.write("#blocks in " + name + ": " + str(len(blocks)) + "\n")

    dists = [[min(seq_identity(str(mn.seq), bl), seq_identity(str(mn.seq), rc(bl))) for mn in mns] for bl in blocks]

    sq_sum = get_sq_sum(blocks, dists)
    f.write("Sq mean for blocks ("  + name + "): "  + str(sq_sum) + "\n")

    global res_str
    res_str += "\t" + "{:.4f}".format(sq_sum)
    return sq_sum


def DBindex_inner(blocks, distfun, mns):
    dists = [[min(distfun(str(mn.seq), bl), distfun(str(mn.seq), rc(bl))) for mn in mns] for bl in blocks]
    blmn = [dists[i].index(min(dists[i])) for i in range(len(blocks))]
    radius = [0] * len(mns)
    for i in range(len(blocks)):
        radius[blmn[i]] = max(radius[blmn[i]], dists[i][blmn[i]])

    DBindex = 0
    for i in range(len(mns)):
        dbmx = 0
        for j in range(len(mns)):
            if i == j:
                continue
            if distfun(mns[i].seq, mns[j].seq) == 0:
                dbmx = 100
                continue
            if distfun(mns[i].seq, rc(mns[j].seq)) == 0:
                dbmx = 100
                continue

            dbmx = max(dbmx,
                       (radius[i] + radius[j]) / min(distfun(mns[i].seq, mns[j].seq),
                                                     distfun(mns[i].seq, rc(mns[j].seq))))

        DBindex += dbmx
    DBindex /= len(mns)
    return DBindex


def DaviesBouldinIndex(path_seq, tsv, mns, f, name):
    blocks = get_blocks(path_seq, tsv)
    f.write("#blocks in " + name + ": " + str(len(blocks)) + "\n")
    DBindex = DBindex_inner(blocks, seq_identity, mns)
    f.write("DBindex for blocks (" + name + "): " + str(DBindex) + "\n")
    f.flush()
    return DBindex

def calc_DBindex_subseq(seq, tsv_B_res, tsv_Ivan_res, CAmn, IAmn, f):
    blocks = get_blocks_interseq(seq, tsv_B_res, tsv_Ivan_res)
    DBindMN = DBindex_inner(blocks, local_dist, CAmn)
    DBindIV = DBindex_inner(blocks, local_dist, IAmn)

    return DBindMN, DBindIV

def main():
    args = parse_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    f = open(os.path.join(args.outdir, "summary"), "w")
    fa = open(os.path.join(args.outdir, "..", "MonomerComparasio.csv"), "a")

    CA_mon = unique(load_fasta(args.camon))
    IVAN_mon = unique([x for x in load_fasta(args.ivanmon) if "H1" in x.id])

    f.write("CA monomers#: " + str(len(CA_mon)) + "\n")
    f.write("IA monomers#: " + str(len(IVAN_mon)) + "\n")

    tsv_res = run_SD(args.camon, args.seq, os.path.join(args.outdir, "CAmon"))
    MonomersNew = ExtractValuableMonomers.getValuableMonomers(args.seq, tsv_res, CA_mon, args.outdir)
    print("Valuable mon cnt:", len(MonomersNew))
    print(MonomersNew)
    SDutils.savemn(os.path.join(args.outdir, "valm.fa"), MonomersNew)
    print(MonomersNew[0].seq)

    MonomersNew_ = getSameMons(MonomersNew, IVAN_mon)
    print(MonomersNew_)
    f.write("MonomerNew_ matched subset#: " + str(len(MonomersNew_)) + "\n")
    MonomersNew_ = addMostFreq(MonomersNew, MonomersNew_, tsv_res, len(IVAN_mon) - len(MonomersNew_))

    print(len(MonomersNew_))
    print(len(IVAN_mon))

    tsv_MN_res = run_SD_mn(MonomersNew_, args.seq, os.path.join(args.outdir, "MN_mon"))
    tsv_Ivan_res = run_SD_mn(IVAN_mon, args.seq, os.path.join(args.outdir, "IVAN_mon"))

    if len(IVAN_mon) > len(MonomersNew_):
        Ivan_mon_ = getSameMons(IVAN_mon, MonomersNew_)
        Ivan_mon_ = addMostFreq(IVAN_mon, Ivan_mon_, tsv_Ivan_res, len(MonomersNew_) - len(Ivan_mon_))
        IVAN_mon = Ivan_mon_
        tsv_Ivan_res = run_SD_mn(IVAN_mon, args.seq, os.path.join(args.outdir, "IVAN_mon_"))

    sqdMN, sqdIV = calc_mean_for_subseq(args.seq, tsv_MN_res, tsv_Ivan_res, MonomersNew_, IVAN_mon, f)
    dbiMN, dbiIV = calc_DBindex_subseq(args.seq, tsv_MN_res, tsv_Ivan_res, MonomersNew_, IVAN_mon, f)

    fa.write(str(len(IVAN_mon)) + "/" + str(len(MonomersNew_)) + "\t")
    fa.write("{:.2f}".format(sqdIV) + "/" + "{:.2f}".format(sqdMN) + "\t")
    fa.write("{:.2f}".format(dbiIV) + "/" + "{:.2f}".format(dbiMN) + "\t")
    fa.write(str(len(CA_mon)) + "\t")
    fa.write(str(len(MonomersNew)) + "\t")
    fa.write("{:.2f}".format(blocks_sqmean(args.seq, tsv_res, MonomersNew, f, "MonomersNew")) + "\t")
    fa.write("{:.2f}".format(DaviesBouldinIndex(args.seq, tsv_res, MonomersNew, f, "MonomersNew")) + "\n")
    fa.close()

    #f.write("\n\nElbow from zero elems\n")
    #elbow_calc(args.seq, tsv_B_res, [], CA_mon, f)

    #f.write("\n\nElbow from B set\n")
    #elbow_calc(args.seq, tsv_B_res, B, CA_mon, f)


    global res_str
    f.write(res_str)
    f.close()

if __name__ == "__main__":
    main()