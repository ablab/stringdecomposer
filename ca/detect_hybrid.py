#!/usr/bin/env python

import os
import sys
import csv
import argparse
import shutil
import edlib

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO

FreqCeiling = 40*23
TotalMonomerBlocks = 0

#Log class, use it, not print
class Log:
    text = ""

    def log(self, s):
        self.text += s + "\n"
        print(s)

    def warn(self, s):
        msg = "WARNING: " + s
        self.text += msg + "\n"
        sys.stdout.write(msg)
        sys.stdout.flush()

    def err(self, s):
        msg = "ERROR: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()
is_RC_alignment = False

def init_monomers(monomers_path):
    #save info about monomers
    monomers_list = []
    for record in SeqIO.parse(monomers_path, "fasta"):
        monomers_list.append(record)
    return monomers_list


def get_record_by_id(mid, monomer_list):
    for record in monomer_list:
        if record.id == mid:
            return record
    return -1


def parse_args():
    parser = argparse.ArgumentParser(description='Shift monomers')
    parser.add_argument('--finalDec', dest="finalDec",
                        help='path to final_decomposition.tsv file from StringDecompoer output')
    parser.add_argument('--monomers', '-mn', dest="monomers",
                        help='path to fasta file with initial monomers')
    parser.add_argument('--outdir', '-o', dest="outdir", help="path to output directory")
    return parser.parse_args()

def rc(seq):
    res_seq = ""
    for i in range(len(seq)):
        pos = len(seq) - i - 1
        if seq[pos] == 'A':
            res_seq += "T"
        elif seq[pos] == 'T':
            res_seq += "A"
        elif seq[pos] == 'C':
            res_seq += "G"
        elif seq[pos] == "G":
            res_seq += "C"
        else:
            res_seq += seq[pos]
    return res_seq


def seq_identity(a, b):
    result = edlib.align(a, b, mode="NW", task="locations")
    if result["editDistance"] == -1:
        return 10**9
    return result["editDistance"] * 100 / max(len(a), len(b))


def is_same(mn1, mn2, mxDiv):
    return seq_identity(str(mn1.seq), str(mn2.seq)) < mxDiv/2


def calc_mn_order_stat(args):
    global TotalMonomerBlocks

    cnt_mn = {}
    db_cnt = {}
    trp_cnt = {}
    RC_cnt = 0
    FR_cnt = 0
    with open(args.finalDec, "r") as f:
        csv_reader = csv.reader(f, delimiter='\t')
        pprev_row = []
        prev_row = []
        for row in csv_reader:
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
                    RC_cnt += 1
                else:
                    FR_cnt += 1

                if row[-1] != '?':
                    TotalMonomerBlocks += 1
                    if mon not in cnt_mn:
                        cnt_mn[mon] = 0
                    cnt_mn[mon] += 1

                if row[-1] != '?' and prev_row[-1] != '?':
                    if (pmon, mon) not in db_cnt:
                        db_cnt[(pmon, mon)] = 0
                    db_cnt[(pmon, mon)] += 1
                    if pprev_row != []:
                        ppmon = pprev_row[1]
                        if ppmon[-1] == "'":
                            ppmon = ppmon[:-1]
                        ppiden = float(pprev_row[4])
                        if pprev_row[-1] != '?':
                            if (ppmon, pmon, mon) not in trp_cnt:
                                trp_cnt[(ppmon, pmon, mon)] = 0
                            trp_cnt[(ppmon, pmon, mon)] += 1

            pprev_row = prev_row
            prev_row = row
    if RC_cnt > FR_cnt:
        global is_RC_alignment
        is_RC_alignment = True
    return cnt_mn, db_cnt


def save_mn(args, shifted_mn):
    outf = os.path.join(args.outdir, "hybrid_mn.fasta")
    with open(outf, "w") as fa:
        for record in shifted_mn:
            SeqIO.write(record, fa, "fasta")


def get_hybrid_len(main_mn, mn1, mn2):
    resDiv = 500
    mn_identity = 100
    bst_res = (0, 0)
    for prfx in range(30, len(mn1.seq)):
        suffix = len(str(main_mn.seq)) - prfx
        if suffix < 30:
            break
        hbr = str(mn1.seq)[:prfx] + str(mn2.seq)[len(mn2.seq) - suffix:]
        cur_identity = seq_identity(hbr, str(main_mn.seq))
        if mn_identity > cur_identity:
            mn_identity = cur_identity
            bst_res =  (prfx, suffix)

    if (mn_identity * 2 <= resDiv):
        return (bst_res[0], bst_res[1], mn_identity)
    return (0, 0, 100)


def get_align(mn1, mn2):
    result = edlib.align(mn1, mn2, mode="NW", task="path")
    return edlib.getNiceAlignment(result, mn1, mn2)


def detect_hybrid_mn(monomers_list, sft_db_cnt, cnt_mn):
    global TotalMonomerBlocks
    global FreqCeiling

    for i in range(len(monomers_list)):
        log.log("hybrid detection for monomer#" + str(i) + " out of " + str(len(monomers_list)))
        bst_hyber_score = 100
        bst_hyber = (0, 0)

        for j in range(len(monomers_list)):
            for g in range(len(monomers_list)):
                if i == j or i == g:
                    continue

                if cnt_mn[monomers_list[j].id] < TotalMonomerBlocks/FreqCeiling:
                    continue

                if cnt_mn[monomers_list[g].id] < TotalMonomerBlocks/FreqCeiling:
                    continue

                if (cnt_mn[monomers_list[j].id] < cnt_mn[monomers_list[i].id]):
                    continue

                if (cnt_mn[monomers_list[g].id] < cnt_mn[monomers_list[i].id]):
                    continue

                hybrid_res = get_hybrid_len(monomers_list[i], monomers_list[j], monomers_list[g])
                if (hybrid_res[0] != 0 and hybrid_res[2] < bst_hyber_score):
                    bst_hyber_score = hybrid_res[2]
                    bst_hyber = (j, g)

        if bst_hyber_score < 100:
            print(bst_hyber, bst_hyber_score)
            cnt1, cnt2, scr = get_hybrid_len(monomers_list[i], monomers_list[bst_hyber[0]], monomers_list[bst_hyber[1]])
            mn1 = monomers_list[bst_hyber[0]].id.split('_')[1]
            mn2 = monomers_list[bst_hyber[1]].id.split('_')[1]

            old_name = monomers_list[i].id
            print(old_name, mn1, cnt1, mn2, cnt2, scr)
            if scr > 5:
                continue
            concatmn = monomers_list[bst_hyber[0]].seq[:cnt1] + monomers_list[bst_hyber[1]].seq[-cnt2:]
            print(concatmn)
            alignCon = get_align(monomers_list[i].seq, concatmn)
            print(alignCon["query_aligned"])
            print(alignCon["matched_aligned"])
            print(alignCon["target_aligned"])
            print()

            align1 = get_align(monomers_list[i], monomers_list[bst_hyber[0]])
            align2 = get_align(monomers_list[i], monomers_list[bst_hyber[1]])

            print(align1["query_aligned"].seq)
            print(align1["matched_aligned"])
            print(align1["target_aligned"].seq)

            print(align2["query_aligned"].seq)
            print(align2["matched_aligned"])
            print(align2["target_aligned"].seq)

            if (mn1 != mn2):
                monomers_list[i].id += "_hybrid_" + mn1 + "(" + str(cnt1) + ")_" + mn2 + "(" + str(cnt2) + ")"
            else:
                monomers_list[i].id += "_variant_" + mn1

            #monomers_list[i].id += "_hybrid_" + mn1 + "_" + mn2
            cnt_mn[monomers_list[i].id] = cnt_mn[old_name]
            kkeys_list = list(sft_db_cnt.keys())
            for key in kkeys_list:
                if key[0] == old_name:
                    sft_db_cnt[(monomers_list[i].id, key[1])] = sft_db_cnt[key]
                if key[1] == old_name:
                    sft_db_cnt[(key[0], monomers_list[i].id)] = sft_db_cnt[key]

    return monomers_list, sft_db_cnt


def main():
    log.log("Start Shift Monomers")
    args = parse_args()
    monomers_list = init_monomers(args.monomers)
    cnt_mn, db_cnt = calc_mn_order_stat(args)
    monomers_list, sft_db_cnt = detect_hybrid_mn(monomers_list, db_cnt, cnt_mn)
    save_mn(args, monomers_list)

if __name__ == "__main__":
    main()
