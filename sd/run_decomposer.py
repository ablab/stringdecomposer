#!/usr/bin/env python

import scripts.common.identities_count_utils as icu
import scripts.common.files_utils as futils

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

import subprocess

import numpy as np; np.random.seed(0)
import pandas as pd

import re
import edlib

import joblib

p = os.path.abspath(__file__)
logreg_file = os.path.join(os.path.dirname(p),
                           'models',
                           'new_ont_logreg_model.sav')
clf = joblib.load(logreg_file)


def run(sequences, monomers, num_threads, scoring, batch_size, raw_file):
    ins, dels, mm, match = scoring.split(",")
    p = os.path.abspath(__file__)
    sd_exec_file = p[:-len("run_decomposer.py")] + "/bin/dp"
    print("Run", sd_exec_file, " with parameters ", sequences, monomers, num_threads, batch_size, scoring, file=sys.stderr)
    with open(raw_file, 'w') as f:
        subprocess.run([sd_exec_file, sequences, monomers, num_threads, batch_size, ins, dels, mm, match], stdout = f, check = True)
    with open(raw_file, 'r') as f:
        raw_decomposition = "".join(f.readlines())
    return raw_decomposition


def main():
    parser = argparse.ArgumentParser(description='Decomposes string into blocks alphabet')
    parser.add_argument('sequences', help='fasta-file with long reads or genomic sequences')
    parser.add_argument('monomers', help='fasta-file with monomers')
    parser.add_argument('-t', '--threads',  help='number of threads (by default 1)', default="1", required=False)
    parser.add_argument('-o', '--out-file',  help='output tsv-file (by default final_decomposition.tsv)', default="./final_decomposition.tsv", required=False)
    parser.add_argument('-i', '--min-identity',  \
                         help='only monomer alignments with percent identity >= MIN_IDENTITY are printed (by default MIN_IDENTITY=0)', type=int, default=0, required=False)
    parser.add_argument('-s', '--scoring', \
                         help='set scoring scheme for SD in the format "insertion,deletion,mismatch,match" (by default "-1,-1,-1,1")', default="-1,-1,-1,1", required=False)
    parser.add_argument('-b', '--batch-size',  help='set size of the batch in parallelization (by default 5000)', type=str, default="5000", required=False)
    parser.add_argument('--fast',  help='doesn\'t generate second best monomer and homopolymer scores', action="store_true")

    args = parser.parse_args()
    raw_decomposition = run(args.sequences, args.monomers, args.threads, args.scoring, args.batch_size, args.out_file[:-len(".tsv")] + "_raw.tsv")
    print("Saved raw decomposition to " + args.out_file[:-len(".tsv")] + "_raw.tsv", file=sys.stderr)

    reads = futils.load_fasta(args.sequences, "map")
    monomers = futils.load_fasta(args.monomers)
    monomers = icu.add_rc_monomers(monomers)
    print("Transforming raw alignments...", file=sys.stderr)
    icu.convert_tsv(raw_decomposition, reads, monomers, args.out_file, clf, int(args.min_identity), args.fast)
    print("Transformation finished. Results can be found in " + args.out_file, file=sys.stderr)


if __name__ == "__main__":
    main()
