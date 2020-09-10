#!/usr/bin/env python

import common.identities_count_utils as icu
import common.files_utils as futils

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

import edlib

import joblib

p = os.path.abspath(__file__)
logreg_file = os.path.join(os.path.dirname(p), "..",
                           'models',
                           'new_ont_logreg_model.sav')
clf = joblib.load(logreg_file)

def add_rc_monomers(monomers):
    res = []
    for m in monomers:
        res.append(m)
        res.append(futils.make_record(m.seq.reverse_complement(), m.name + "'", m.id + "'"))
    return res


def print_read(fout, fout_alt, dec, read, monomers):
    dec = icu.convert_read(dec, read, monomers, clf)
    for d in dec:
        fout.write("\t".join([read.name, d["m"], d["start"], d["end"], "{:.2f}".format(d["score"]), \
                                                d["second_best"], "{:.2f}".format(d["second_best_score"]), d["q"]]) + "\n")
        for a in d["alt"]:
            star = "-"
            if a == d["m"]:
                star = "*"
            fout_alt.write("\t".join([read.name, a, d["start"], d["end"], "{:.2f}".format(d["alt"][a]), star]) + "\n")


def convert_tsv(filename, reads, monomers, outfile):
    with open(outfile[:-len(".tsv")] + "_alt.tsv", "w") as fout_alt:
        with open(outfile, "w") as fout:
            with open(filename, "r") as fin:
                cur_dec = []
                prev_read = None
                for ln in fin.readlines():
                    read, monomer, start, end = ln.split("\t")[:4]
                    if read != prev_read and prev_read != None:
                        print_read(fout, fout_alt, cur_dec, reads[prev_read], monomers)
                        cur_dec = []
                    prev_read = read
                    start, end = int(start), int(end)
                    cur_dec.append({"m": monomer, "start": start, "end": end})
                if len(cur_dec) > 0:
                    print_read(fout, fout_alt, cur_dec, reads[prev_read], monomers)
                    
def convert_fasta(filename, reads, monomers, outfile):
    with open(outfile[:-len(".tsv")] + "_alt.tsv", "w") as fout_alt:
        with open(outfile, "w") as fout:
            with open(filename, "r") as fin:
                cur_dec = []
                prev_read = None
                for ln in fin.readlines():
                    if ln.startswith(">"):
                        read = ln.split("/")[0][1:]
                        if read != prev_read and prev_read != None:
                            print_read(fout, fout_alt, cur_dec, reads[prev_read], monomers)
                            cur_dec = []
                        prev_read = read
                        start, end = [int(x) for x in ln.split("/")[1].split("_")]
                        cur_dec.append({"m": None, "start": start, "end": end})
                if len(cur_dec) > 0:
                    print_read(fout, fout_alt, cur_dec, reads[prev_read], monomers)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert decomposition scores into identities')
    parser.add_argument('-s', '--sequences', help='fasta-file with long reads sequences', required=True)
    parser.add_argument('-m', '--monomers', help='fasta-file with monomers', required=True)
    parser.add_argument('-d', '--decomposition', help='tsv-file (for DP) or fasta-file (for AC) with decomposition', required=True)
    parser.add_argument('-o', '--out',  help='output tsv-file, by default will be saved into decomposition.tsv', required=False)

    args = parser.parse_args()
    outfile = args.out
    if outfile == None:
        outfile = "./decomposition.tsv"

    reads = futils.load_fasta(args.sequences, "map")
    monomers = futils.load_fasta(args.monomers)
    monomers = add_rc_monomers(monomers)
    if args.decomposition.endswith("tsv"):
        convert_tsv(args.decomposition, reads, monomers, outfile)
    else:
        convert_fasta(args.decomposition, reads, monomers, outfile)