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


def parse_args():
    parser = argparse.ArgumentParser(description='Compare Origin and Reported monomers')
    parser.add_argument('MGdir', help='path to dir with MonomerGenerator results')
    parser.add_argument('monomers', help='fasta-file with origin monomers')
    return parser.parse_args()


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
                a = str(origin_monomers[omon].seq)
                b = str(reported_monomers[emon].seq)
                result = edlib.align(a, b, mode="NW", task="locations")
                dist = result["editDistance"] * 100 / max(len(a), len(b))

                result = edlib.align(a, b, mode="NW", task="path")
                nice = edlib.getNiceAlignment(result, a, b)
                writer.writerow([omon, emon, str(dist), "\n".join(nice.values())])

if __name__ == "__main__":
    main()