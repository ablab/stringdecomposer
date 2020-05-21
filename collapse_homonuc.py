#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import os
from os import listdir
from os.path import join, isfile
import sys
import argparse

from itertools import groupby


def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
        for r in records:
            records[r] = records[r].upper() 
    else:
        records = list(SeqIO.parse(filename, "fasta"))
        for i in range(len(records)):
            records[i] = records[i].upper()
    return records

def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)

def save_fasta(filename, orfs):
    with open(filename, "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")

def convert_to_homo(seq):
    res = ''.join(i for i, _ in groupby(seq))
    return res

if __name__ == "__main__":
	filename = sys.argv[1]
	reads = load_fasta(filename)
	res = [make_record(Seq(convert_to_homo(x.seq)), x.name, x.id) for x in reads]
	save_fasta(".".join(filename.split(".")[:-1]) + "_homocollapsed.fa", res)
	print("Saved to", ".".join(filename.split(".")[:-1]) + "_homocollapsed.fa")
