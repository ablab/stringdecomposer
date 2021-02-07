#!/usr/bin/env python

import os
import sys
import csv
import argparse
import shutil
import edlib

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

FreqCeiling = 40
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
        record.id = "_".join(record.id.split("_")[:2])
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
    parser.add_argument('--shift', dest="shift", type=int, help='shift for monomers (Monomer << shift)')
    parser.add_argument("--RC", dest="RC", action='store_true')
    parser.add_argument('--outdir', '-o', dest="outdir", help="path to output directory")
    return parser.parse_args()


def get_edge_list(args):
    mongr_path = os.path.join(args.MGres, "final", "monomer_gr")

    edge_list = []
    with open(mongr_path, 'r') as f:
        for line in f:
            v1, v2, wg = line.split()
            wg = int(wg)
            edge_list.append((v1, v2, wg))
    return edge_list


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


def get_one_shifted_monomer(mn1, mn2, shift, args):
    global is_RC_alignment
    edge_name = "mn_" + mn1.id.split('_')[-1] + "-" + mn2.id.split('_')[-1]
    seq1 = mn1.seq
    seq2 = mn2.seq
    if (args.RC):
        seq1 = rc(str(mn1.seq))
        seq2 = rc(str(mn2.seq))

    nseq = str(seq1)[shift:] + str(seq2)[:shift]
    if is_RC_alignment != args.RC:
        nseq = str(seq2)[shift:] + str(seq1)[:shift]

    mn_record = SeqRecord(Seq(nseq), id=edge_name, description="")
    return mn_record


def generate_shifted_monomers(monomers_list, shift, db_cnt, trp_cnt, args):
    log.log("Start generate shifted monomers")
    sft_db_cnt = {}
    sh_monomers = []
    cnt_mn = {}
    for edge in db_cnt:
        if edge in db_cnt and db_cnt[edge] > 0:
            mn1 = get_record_by_id(edge[0], monomers_list)
            mn2 = get_record_by_id(edge[1], monomers_list)
            sh_monomers.append((get_one_shifted_monomer(mn1, mn2, shift, args), edge[0], edge[1]))
            cnt_mn[sh_monomers[-1][0].id] = db_cnt[edge]

    for shmn1 in sh_monomers:
        for shmn2 in sh_monomers:
            if shmn1[2] == shmn2[1]:
                if (shmn1[0].id, shmn2[0].id) not in sft_db_cnt:
                    sft_db_cnt[(shmn1[0].id, shmn2[0].id)] = 0
                sft_db_cnt[(shmn1[0].id, shmn2[0].id)] = 0
                if (shmn1[1], shmn1[2], shmn2[2]) in trp_cnt:
                    sft_db_cnt[(shmn1[0].id, shmn2[0].id)] = trp_cnt[(shmn1[1], shmn1[2], shmn2[2])]
    res_sh_mn = []
    for elem in sh_monomers:
        res_sh_mn.append(elem[0])
    return res_sh_mn, sft_db_cnt, cnt_mn


def seq_identity(a, b):
    result = edlib.align(a, b, mode="NW", task="locations")
    if result["editDistance"] == -1:
        return 10**9
    return result["editDistance"] * 100 / max(len(a), len(b))


def is_same(mn1, mn2, mxDiv):
    return seq_identity(str(mn1.seq), str(mn2.seq)) < mxDiv/2


def delete_same_mn(shifted_mn, sft_db_cnt, cnt_mn):
    log.log("Start delete same monomers")
    mxDiv = 5
    delete = [0] * len(shifted_mn)
    for i in range(len(shifted_mn)):
        print(i, delete)
        mn1 = shifted_mn[i]
        for j in range(i + 1, len(shifted_mn)):
            if delete[i] == 1 or delete[j] == 1:
                continue

            mn2 = shifted_mn[j]
            if is_same(mn1, mn2, mxDiv):
                if cnt_mn[mn1.id] > cnt_mn[mn2.id]:
                    delete[j] = 1
                else:
                    delete[i] = 1
                    mn1, mn2 = mn2, mn1

                for g in range(len(shifted_mn)):
                    if delete[g] == 0:
                        mn3 = shifted_mn[g]
                        if (mn1.id, mn3.id) not in sft_db_cnt:
                            sft_db_cnt[(mn1.id, mn3.id)] = 0
                        if (mn2.id, mn3.id) not in sft_db_cnt:
                            sft_db_cnt[(mn2.id, mn3.id)] = 0
                        if (mn3.id, mn1.id) not in sft_db_cnt:
                            sft_db_cnt[(mn3.id, mn1.id)] = 0
                        if (mn3.id, mn2.id) not in sft_db_cnt:
                            sft_db_cnt[(mn3.id, mn2.id)] = 0

                        sft_db_cnt[(mn1.id, mn3.id)] += sft_db_cnt[(mn2.id, mn3.id)]
                        sft_db_cnt[(mn3.id, mn1.id)] += sft_db_cnt[(mn3.id, mn2.id)]
    res_sh_mn = []
    for i in range(len(shifted_mn)):
        if delete[i] == 0:
            res_sh_mn.append(shifted_mn[i])

    return res_sh_mn, sft_db_cnt


def calc_mn_order_stat(args):
    global TotalMonomerBlocks

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
    return db_cnt, trp_cnt


def save_init_graph(args, db_cnt, monomers_list, outfile):
    dotst = "digraph MonomerGraph {\n"

    for i in range(len(monomers_list)):
        dotst += "\"" + monomers_list[i].id + "\";\n"

    for i in range(len(monomers_list)):
        curm = monomers_list[i].id
        for j in range(len(monomers_list)):
            m2 = monomers_list[j].id
            if (curm, m2) in db_cnt and db_cnt[(curm, m2)] > 0:
                dotst += "\"" + curm + "\" -> \"" + m2 + "\" [label=\"" + str(db_cnt[(curm, m2)]) + "\"];\n"


    dotst += "}\n"

    dotpath = os.path.join(args.outdir, outfile)
    with open(dotpath, "w") as fw:
        fw.write(dotst)


def save_mn(args, shifted_mn):
    outf = os.path.join(args.outdir, "shifted_mn.fasta")
    with open(outf, "w") as fa:
        for record in shifted_mn:
            SeqIO.write(record, fa, "fasta")


def get_hybrid_len(main_mn, mn1, mn2):
    resDiv = 5
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
            print(old_name, mn1, cnt1, mn2, cnt2)

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
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    monomers_list = init_monomers(args.monomers)
    db_cnt, trp_cnt = calc_mn_order_stat(args)
    save_init_graph(args, db_cnt, monomers_list, "init_monomer_graph.dot")
    shifted_mn, sft_db_cnt, cnt_mn = generate_shifted_monomers(monomers_list, args.shift, db_cnt, trp_cnt, args)
    shifted_mn, sft_db_cnt = delete_same_mn(shifted_mn, sft_db_cnt, cnt_mn)
    shifted_mn, sft_db_cnt = detect_hybrid_mn(shifted_mn, sft_db_cnt, cnt_mn)
    save_init_graph(args, sft_db_cnt, shifted_mn, "final_monomer_graph.dot")
    save_mn(args, shifted_mn)


if __name__ == "__main__":
    main()
