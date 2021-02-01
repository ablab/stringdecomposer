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

class MonomericBlock(object):
    read_name = ""
    lft = 0
    rgh = 0
    seq = ""

    def __init__(self, read_name="", lft=0, rgh=0, seq=""):
        self.read_name = read_name
        self.lft = lft
        self.rgh = rgh
        self.seq = seq


def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
        for r in records:
            #print("Read_name in map: ", r)
            records[r] = records[r].upper()
    else:
        records = list(SeqIO.parse(filename, "fasta"))
        for i in range(len(records)):
            records[i] = records[i].upper()
    return records


def load_last_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    records[-1] = records[-1].upper()
    return records[-1]


def get_min_dist(last_monomer, origin_monomers):
    res_dist = 100
    for omon in origin_monomers:
        a = str(origin_monomers[omon].seq)
        b = str(last_monomer.seq)
        result = edlib.align(a, b, mode="NW", task="locations")
        dist = result["editDistance"]# * 100 / max(len(a), len(b))
        res_dist = min(res_dist, dist)

    return res_dist


def parse_args():
    parser = argparse.ArgumentParser(description='Compare Origin and Reported monomers')
    parser.add_argument('MGdir', help='path to dir with MonomerGenerator results')
    parser.add_argument('monomers', help='fasta-file with origin monomers')
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


def set_blocks_seq(sequences, blocks):
    records = load_fasta(sequences, tp="map")
    for i in range(len(blocks)):
        blocks[i].seq = records[blocks[i].read_name][blocks[i].lft:blocks[i].rgh + 1]


def getMonomerBlocks(mnid, res_tsv, sequences):
    resDiv = 5
    monomer_resolved = []

    with open(res_tsv, "r") as f:
        csv_reader = csv.reader(f, delimiter='\t')
        for row in csv_reader:
            if row[2] == "start":
                continue
            identity = float(row[4])
            if row[-1] == '?':
                continue

            #print("Identity: ", identity)
            #if identity >= 100 - resDiv:
            if row[1][-1] == "'":
                row[1] = row[1][:-1]
                continue

            if (row[1] == mnid):
                monomer_resolved.append(MonomericBlock(row[0], int(row[2]), int(row[3])))

    set_blocks_seq(sequences, monomer_resolved)
    return monomer_resolved


def get_radius(last_monomer, monomers_blocks):
    dists_vector = []
    max_dist = 0
    for mnblock in monomers_blocks:
        a = str(mnblock.seq.seq)
        b = str(last_monomer.seq)
        result = edlib.align(a, b, mode="NW", task="locations")
        dists_vector.append(result["editDistance"])
        max_dist = max(max_dist, result["editDistance"])
    print("Monomers cnt:", len(monomers_blocks))
    dists_vector.sort()
    print("Dists sorted:", dists_vector)
    return  max_dist


def calc_radius(last_monomer, final_dec_path, seq):
    monomers_blocks = getMonomerBlocks(last_monomer.id,
                                       final_dec_path,
                                       seq)
    return get_radius(last_monomer, monomers_blocks)


def resolved_cnt(monomer, final_dec_path, sequences):
    return len(getMonomerBlocks(monomer, final_dec_path, sequences))


def cal_separation(monomer, reported_monomers, final_dec_path, sequences, thr=500):
    res_dist = 100
    sepName = ""
    for omon in reported_monomers:
        if reported_monomers[omon].id == monomer.id:
            continue
        #if resolved_cnt(reported_monomers[omon].id, final_dec_path, sequences) < thr:
        #    continue

        a = str(reported_monomers[omon].seq)
        b = str(monomer.seq)
        result = edlib.align(a, b, mode="NW", task="locations")
        dist = result["editDistance"]# * 100 / max(len(a), len(b))
        res_dist = min(res_dist, dist)
        if res_dist == dist:
            sepName = omon

        b = rc(str(monomer.seq))
        result = edlib.align(a, b, mode="NW", task="locations")
        dist = result["editDistance"]  # * 100 / max(len(a), len(b))
        res_dist = min(res_dist, dist)
        if res_dist == dist:
            sepName = omon

    return res_dist, sepName


def get_cnt_in_centromers(res_tsv, mnid):
    resDiv = 5
    cnt_cet = {}

    with open(res_tsv, "r") as f:
        csv_reader = csv.reader(f, delimiter='\t')
        for row in csv_reader:
            if row[2] == "start":
                continue
            identity = float(row[4])
            if row[-1] == '?':
                continue

            #print("Identity: ", identity)
            if identity < 100 - resDiv:
                continue

            if row[1][-1] == "'":
                row[1] = row[1][:-1]
                continue

            if (row[1] == mnid):
                if row[0].split(':')[0] not in cnt_cet:
                    cnt_cet[row[0].split(':')[0]] = 1
                cnt_cet[row[0].split(':')[0]] += 1
    res_lst = []
    for key in cnt_cet:
        res_lst.append((-cnt_cet[key], key))
    res_lst.sort()
    res_cen, res_cnt = [], []
    for cnt1, cen1 in res_lst:
        res_cen.append(cen1)
        res_cnt.append(-cnt1)
    return res_cen, res_cnt


def main():
    args = parse_args()
    origin_monomers = load_fasta(args.monomers, "map")
    reported_monomers = load_fasta(os.path.join(args.MGdir, "final", "monomers.fa"), "map")
    frq_monomers = {}
    cmpMonomers = os.path.join(args.MGdir, "cmpReportedAndOriginMonomers.csv")
    with open(cmpMonomers, "w") as fw:
         writer = csv.writer(fw)
         writer.writerow(['OriginMonomerName', 'ReportedMonomerName', 'Distance', 'Radius', 'Separation', 'Alignment'])
         for emon in reported_monomers:
             radius = 0#calc_radius(reported_monomers[emon],
                        #          os.path.join(args.MGdir, "final", "final_decomposition.tsv"),
                        #          os.path.join(args.MGdir, "sequence.fa"))
             separation = 0#cal_separation(reported_monomers[emon], reported_monomers, os.path.join(args.MGdir, "final", "final_decomposition.tsv"),
                            #     os.path.join(args.MGdir, "sequence.fa"))

             for omon in origin_monomers:
                 a = str(origin_monomers[omon].seq)
                 b = str(reported_monomers[emon].seq)

                 result = edlib.align(a, b, mode="NW", task="locations")
                 dist1 = result["editDistance"] * 100 / max(len(a), len(b))
                 result = edlib.align(rc(a), b, mode="NW", task="locations")
                 dist2 = result["editDistance"] * 100 / max(len(a), len(b))
                 if dist2 < dist1:
                     a = rc(a)

                 result = edlib.align(a, b, mode="NW", task="locations")
                 dist = result["editDistance"] * 100 / max(len(a), len(b))

                 result = edlib.align(a, b, mode="NW", task="path")
                 nice = edlib.getNiceAlignment(result, a, b)
                 writer.writerow([omon, emon, str(dist), str(radius), str(separation), "\n".join(nice.values())])

    summary_path = os.path.join(args.MGdir, "summary.csv")
    summary_upgrade = os.path.join(args.MGdir, "summary_upgrade.csv")

    with open(summary_path, "r") as f:
        reader = csv.reader(f)
        for line in reader:
            if line[4][0] == 'S':
                continue
            if int(line[4]) >= 400:
                frq_monomers["mn_" + line[0]] = reported_monomers["mn_" + line[0]]

    with open(summary_path, "r") as f:
        reader = csv.reader(f)
        with open(summary_upgrade, "w") as fw:
            writer = csv.writer(fw)
            line_id = 0
            for line in reader:
                print(line)
                if line_id == 0:
                    writer.writerow(line + ["Dist to original monomers", "Dist to previously generated monomer", "Length", "Centromers position", "Number of occurance", "Dist to all other monomers", "Dist to frequent monomers"])
                else:
                    nxt_folder_mono = os.path.join(args.MGdir, "iter_"  + str(line_id), "monomers.fa")
                    final_monomers_folder = os.path.join(args.MGdir, "final", "monomers.fa")
                    mn_name = "mn_" + str(line_id - 1)
                    if True:#(os.path.isfile(nxt_folder_mono)):
                        last_monomer = load_fasta(final_monomers_folder, "map")[mn_name]
                        nxt_folder_mono = os.path.join(args.MGdir, "iter_" + str(line_id - 1), "monomers.fa")
                        all_monomers = load_fasta(final_monomers_folder, "map")
                        nw_monomers = {}
                        for i in range(line_id - 1):
                            if ("mn_" + str(i)) in all_monomers:
                                nw_monomers["mn_" + str(i)] = all_monomers["mn_" + str(i)]

                        print(nw_monomers)
                        mdist = get_min_dist(last_monomer, origin_monomers)
                        #del nw_monomers[last_monomer.id]
                        mgdist = get_min_dist(last_monomer, nw_monomers)
                        centromers_list, cnt_in_cen = get_cnt_in_centromers(os.path.join(args.MGdir, "final", "final_decomposition.tsv"), mn_name)

                        separation, sepName1 = cal_separation(last_monomer, all_monomers, os.path.join(args.MGdir, "final", "final_decomposition.tsv"),
                                                    os.path.join(args.MGdir, "sequence.fa"), 0)

                        separationFrq, sepName2 = cal_separation(last_monomer, frq_monomers, os.path.join(args.MGdir, "final", "final_decomposition.tsv"),
                                                    os.path.join(args.MGdir, "sequence.fa"), 50)
                        writer.writerow(line + [str(mdist), str(mgdist), str(len(last_monomer.seq)), str(centromers_list), str(cnt_in_cen), str(separation), str(separationFrq), sepName1, sepName2])
                    else:
                        writer.writerow(line + ["-", "-", "-", "-", "-"])
                line_id += 1
                fw.flush()

if __name__ == "__main__":
    main()