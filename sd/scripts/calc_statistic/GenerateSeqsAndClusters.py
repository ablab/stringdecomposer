#/usr/bin/env python3

import argparse
import edlib
import os
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO
import pandas as pd
import numpy as np

mn_names = ["S3C17H1-B.7", "S3C17H1L.7", "S3C17H1-C.7"]
#mn_names = ["mn_43", "mn_11"]


def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-sd_tsv")
    parser.add_argument("-seq")
    parser.add_argument("-seq_o")
    parser.add_argument("-cluster_o")

    return parser.parse_args()

def main():
    args = parse_args()

    seqs_dict = {}
    for record in SeqIO.parse(args.seq, "fasta"):
        seqs_dict[record.id] = str(record.seq).upper()

    df_sd = pd.read_csv(args.sd_tsv, "\t")
    print(df_sd.iloc[:, 4].head())
    df_sd.iloc[:, 4] = df_sd.iloc[:, 4].apply(int)
    df_sd = df_sd.loc[df_sd.iloc[:, 4] > 60]
    df_sd = df_sd.loc[df_sd.iloc[:, 1].isin(mn_names)]

    with open(args.seq_o,"w") as seqf:
        with open(args.cluster_o, "w") as clsf:
            for i in range(len(df_sd)):
                record = SeqRecord(Seq(seqs_dict[df_sd.iloc[i,0]][df_sd.iloc[i,2]:(df_sd.iloc[i,3] + 1)]),
                                   id="bl_" + str(i), description="")
                SeqIO.write(record, seqf, "fasta")
                clsf.write(record.id + " " + str(mn_names.index(df_sd.iloc[i, 1])) + "\n")

if __name__ == "__main__":
    main()