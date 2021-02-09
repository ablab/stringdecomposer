from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import sys
import os
from os import listdir
from os.path import isfile, isdir, join

import re
import edlib
import random

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please specify tsv-file to convert: python3 convert2bed.py <filename>")
        exit(-1)
    tsvfile = sys.argv[1]
    outfile = tsvfile[:-len(".tsv")] + ".bed"
    mono = {}
    with open(outfile, "w") as fout:
        with open(tsvfile, "r") as fin:
            for ln in fin.readlines():
                lst = ln.strip().split("\t")
                tig_name = lst[0].split(":")[0]
                tstart, tend = lst[0].split(":")[1].split("-")
                tstart, tend = int(tstart), int(tend)
                rev = "+"
                if lst[1].endswith("'"):
                    lst[1] = lst[1][:-1]
                    rev = "-"
                if lst[1] in mono:
                    r,g,b = mono[lst[1]]
                else:
                    r,g,b = random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)
                    mono[lst[1]] = [r,g,b]
                fout.write("\t".join([tig_name, str(tstart + int(lst[2])), str(tstart + int(lst[3]))\
                                    , lst[1], str(int(float(lst[4]))), rev, str(tstart + int(lst[2])), str(tstart + int(lst[3]))\
                                    , str(r) + "," + str(g) + "," + str(b)]) + "\n")
    print("Please find result in ", outfile)