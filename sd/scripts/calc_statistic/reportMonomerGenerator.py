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


def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
        for r in records:
            print("Read_name in map: ", r)
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
        dist = result["editDistance"] * 100 / max(len(a), len(b))
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

def main():
    args = parse_args()
    origin_monomers = load_fasta(args.monomers, "map")
    reported_monomers = load_fasta(os.path.join(args.MGdir, "monomers.fa"), "map")
    cmpMonomers = os.path.join(args.MGdir, "cmpReportedAndOriginMonomers.csv")
    with open(cmpMonomers, "w") as fw:
        writer = csv.writer(fw)
        writer.writerow(['OriginMonomerName', 'ReportedMonomerName', 'Distance', 'Alignment'])
        for omon in origin_monomers:
            for emon in reported_monomers:
                a = rc(str(origin_monomers[omon].seq))
                b = rc(str(reported_monomers[emon].seq))
                result = edlib.align(a, b, mode="NW", task="locations")
                dist = result["editDistance"] * 100 / max(len(a), len(b))

                result = edlib.align(a, b, mode="NW", task="path")
                nice = edlib.getNiceAlignment(result, a, b)
                writer.writerow([omon, emon, str(dist), "\n".join(nice.values())])

    summary_path = os.path.join(args.MGdir, "summary.csv")
    summary_upgrade = os.path.join(args.MGdir, "summary_upgrade.csv")
    with open(summary_path, "r") as f:
        reader = csv.reader(f)
        with open(summary_upgrade, "w") as fw:
            writer = csv.writer(fw)
            line_id = 0
            for line in reader:
                if line_id == 0:
                    writer.writerow(line + ["Dist to original monomers"])
                else:
                    nxt_folder_mono = os.path.join(args.MGdir, "iter_"  + str(line_id), "monomers.fa")
                    if (os.path.isfile(nxt_folder_mono)):
                        last_monomer = load_last_fasta(nxt_folder_mono)
                        mdist = get_min_dist(last_monomer, origin_monomers)
                        writer.writerow(line + [str(mdist)])
                    else:
                        writer.writerow(line + ["-"])
                line_id += 1

if __name__ == "__main__":
    main()