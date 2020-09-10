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
    monomers = icu.add_rc_monomers(monomers)
    if args.decomposition.endswith("tsv"):
        with open(args.decomposition, 'r') as f:
            raw_decomposition = "".join(f.readlines())
            icu.convert_tsv(raw_decomposition, reads, monomers, outfile, clf, identity_th=0, light=False)
    else:
        icu.convert_fasta(args.decomposition, reads, monomers, outfile, clf)