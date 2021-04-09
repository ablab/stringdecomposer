#!/usr/bin/env python3

import argparse
import edlib
import os
import shutil
import sys
import pandas as pd
import numpy as np
import SDutils
from SDutils import rc
from SDutils import sys_call
from SDutils import unique
from SDutils import load_fasta
from SDutils import run_SD
import TriplesMatrix

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--cenID", dest="cenID")
    parser.add_argument("--mon", dest="mon")
    parser.add_argument("--seq", dest="seq")
    parser.add_argument("-o", dest = "outdir")

    return parser.parse_args()


def getMnSim(mon):
    sim = {}
    for mn1 in mon:
        for mn2 in mon:
            sim[(mn1.id, mn2.id)] = min(SDutils.seq_identity(mn1.seq, mn2.seq),
                                        SDutils.seq_identity(mn1.seq, rc(mn2.seq)))
    return sim


cntMerge = 0
cntSplit = 0
tsv_res = ""

def save_seqs(blocks, cluster_seqs_path):
    with open(cluster_seqs_path, "w") as fa:
        for i in range(len(blocks)):
            name = "block" + str(i)
            new_record = SeqRecord(Seq(blocks[i]), id=name, name=name, description="")
            SeqIO.write(new_record, fa, "fasta")

def get_consensus_seq(cluster_seqs_path, arg_threads):
    from Bio.Align.Applications import ClustalwCommandline
    from Bio.Align.Applications import ClustalOmegaCommandline
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Align import MultipleSeqAlignment

    aln_file = '.'.join(cluster_seqs_path.split('.')[:-1]) + "_aln.fasta"
    cmd = ClustalOmegaCommandline(infile=cluster_seqs_path, outfile=aln_file, force=True, threads=arg_threads)
    stdout, stderr = cmd()
    align = AlignIO.read(aln_file, "fasta")

    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.gap_consensus(threshold=0, ambiguous='N')
    consensus = str(consensus).replace('-', '')
    return consensus


def MergeMonomers(mn1, mn2, odir, mons, path_seq):
    print("====MERGE===" + mn1.id + "+" + mn2.id)
    global tsv_res
    resmns = [mn for mn in mons if mn.id != mn1.id and mn.id != mn2.id]

    blocks = SDutils.get_blocks(path_seq, tsv_res, [mn1.id, mn2.id])
    save_seqs(blocks, os.path.join(odir, "blseq.fa"))
    consensus = get_consensus_seq(os.path.join(odir, "blseq.fa"), 16)
    name = mn1.id + "+" + mn2.id
    new_record = SeqRecord(Seq(consensus), id=name, name=name, description="")
    resmns.append(new_record)
    return resmns


def get_blocks(trpl, path_seq, tsv_res):
    block1, block2 = [],[]
    seqs_dict = {}
    for record in SeqIO.parse(path_seq, "fasta"):
        seqs_dict[record.id] = str(record.seq).upper()

    df_sd = pd.read_csv(tsv_res, "\t")
    for i in range(1, len(df_sd) - 1):
        if df_sd.iloc[i,4] > 60:
            if df_sd.iloc[i, 1].rstrip("'") == trpl[1]:
                #print(df_sd.iloc[i - 1, 1], df_sd.iloc[i, 1], df_sd.iloc[i + 1, 1])
                #print(trpl)
                if df_sd.iloc[i + 1, 1].rstrip("'") == trpl[0] and df_sd.iloc[i - 1, 1].rstrip("'") == trpl[2]:
                    block2.append(seqs_dict[df_sd.iloc[i,0]][df_sd.iloc[i,2]:(df_sd.iloc[i,3] + 1)])
                    if df_sd.iloc[i, 1][-1] == "'":
                        block2[-1] = rc(block2[-1])
                else:
                    block1.append(seqs_dict[df_sd.iloc[i, 0]][df_sd.iloc[i, 2]:(df_sd.iloc[i, 3] + 1)])
                    if df_sd.iloc[i, 1][-1] == "'":
                        block1[-1] = rc(block1[-1])
    return block1, block2


def SplitMn(trpl, nnm, odir, mons, path_seq):
    print("====SPLIT===", trpl, nnm)
    global tsv_res
    resmns = [mn for mn in mons if mn.id != trpl[1]]

    block1, block2 = get_blocks(trpl, path_seq, tsv_res)
    save_seqs(block1, os.path.join(odir, "blseq.fa"))
    consensus = get_consensus_seq(os.path.join(odir, "blseq.fa"), 16)
    name = trpl[1]
    new_record = SeqRecord(Seq(consensus), id=name, name=name, description="")
    resmns.append(new_record)

    save_seqs(block2, os.path.join(odir, "blseq.fa"))
    consensus = get_consensus_seq(os.path.join(odir, "blseq.fa"), 16)
    name = nnm
    new_record = SeqRecord(Seq(consensus), id=name, name=name, description="")
    resmns.append(new_record)

    return resmns


def Iteration(iterNum, args, monsPath):
    global tsv_res
    global cntMerge
    global cntSplit

    mons = unique(load_fasta(monsPath))
    sim = getMnSim(mons)

    odir = os.path.join(args.outdir, "i" + str(iterNum))
    if not os.path.exists(odir):
        os.makedirs(odir)

    tsv_res = run_SD(monsPath, args.seq, os.path.join(odir, "InitSD"))
    k_cnt = TriplesMatrix.calc_mn_order_stat(tsv_res, "cen" + str(args.cenID) + "_", maxk=3)
    posscore, cenvec = TriplesMatrix.handleAllMn(k_cnt[1], k_cnt[0], thr=0)

    def get_best(posscore):
        bstm = (-1, -1)
        bp = 0
        for i in range(len(mons)):
            for j in range(len(mons)):
                mn1 = mons[i]
                mn2 = mons[j]
                if i == j:
                    continue
                if sim[(mn1.id, mn2.id)] < 6:
                    scr = posscore[(mn1.id, mn2.id)]
                    if scr > bp:
                        bp = scr
                        bstm = (i, j)
        return bstm, bp

    bstm, bp = get_best(posscore)

    if bp > 0.3:
        cntMerge += 1
        resmn = MergeMonomers(mons[bstm[0]], mons[bstm[1]], odir, mons, args.seq)
        SDutils.savemn(os.path.join(odir, "mn.fa"), resmn)
        return True

    posscore, cenvec = TriplesMatrix.PrefixPosScore(k_cnt[1], k_cnt[0], thr=0)
    bstm, bp = get_best(posscore)
    if bp > 0.3:
        cntMerge += 1
        resmn = MergeMonomers(mons[bstm[0]], mons[bstm[1]], odir, mons, args.seq)
        SDutils.savemn(os.path.join(odir, "mn.fa"), resmn)
        return True

    posscore, cenvec = TriplesMatrix.SuffixPosScore(k_cnt[1], k_cnt[0], thr=0)
    bstm, bp = get_best(posscore)
    if bp > 0.3:
        cntMerge += 1
        resmn = MergeMonomers(mons[bstm[0]], mons[bstm[1]], odir, mons, args.seq)
        SDutils.savemn(os.path.join(odir, "mn.fa"), resmn)
        return True

    exchTrp = TriplesMatrix.SplitAllMn(k_cnt[1], k_cnt[0])
    if len(exchTrp) > 0:
        for key in exchTrp.keys():
            cntSplit += 1
            resmn = SplitMn(key, exchTrp[key], odir, mons, args.seq)
            SDutils.savemn(os.path.join(odir, "mn.fa"), resmn)
            return True
    return False


def main():
    args = parse_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    fa = open(os.path.join(args.outdir, "..", "MergeSplitStat.csv"), "a")
    fa.write(str(args.cenID) + " ")
    fa.write(str(len(unique(load_fasta(args.mon)))) + " ")
    cmonPath = args.mon
    shutil.copyfile(cmonPath, os.path.join(args.outdir, "mn.fa"))

    iterNum = 0
    while(Iteration(iterNum, args, cmonPath)):
        iterNum += 1
        cmonPath = os.path.join(args.outdir, "i" + str(iterNum - 1), "mn.fa")
        shutil.copyfile(cmonPath, os.path.join(args.outdir, "mn.fa"))

    print(os.path.join(args.outdir, "i" + str(iterNum), "InitSD" ,"final_decomposition.tsv"))
    shutil.copyfile(os.path.join(args.outdir, "i" + str(iterNum), "InitSD" ,"final_decomposition.tsv"), os.path.join(args.outdir, "fdec.tsv"))


    fa.write(str(cntMerge) + " " + str(cntSplit) + " ")
    finalm = unique(load_fasta(cmonPath))

    fa.write(str(len(finalm)) + " ")
    fa.write(str(SDutils.blocks_sqmean(args.seq, tsv_res , finalm)) + " ")
    fa.write(str(SDutils.DaviesBouldinIndex(args.seq, tsv_res , finalm)) + "\n")

    fa.close()

if __name__ == "__main__":
    main()