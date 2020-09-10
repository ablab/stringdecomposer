#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord


def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)


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


def save_fasta(filename, orfs):
    with open(filename, "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")