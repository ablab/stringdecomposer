#!/usr/bin/env python

# requirements: clustalw2
import os
import sys
import csv
import argparse
import shutil
import edlib
import multiprocessing
from random import *

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO

from threading import Thread

#Log class, use it, not print
class Log:
    text = ""

    def log(self, s):
        self.text += s + "\n"
        print(s)

    def warn(self, s):
        msg = "WARNING: " + s
        self.text += msg + "\n"
        sys.stdout.write(msg)
        sys.stdout.flush()

    def err(self, s):
        msg = "ERROR: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()


def parse_args():
    parser = argparse.ArgumentParser(description='Monomer Inference Problem: complement monomers set')
    parser.add_argument('-in', '--in-dir', dest="indir", help='directory with monomers and summary (default=\'.\')', default=".",
                        required=False)
    parser.add_argument('-o', '--out-dir', dest="outdir", help='output directory for separate mon (default=\'.\')', default=".", required=False)
    return parser.parse_args()


def init_monomers(monomers_path):
    #save info about monomers
    monomers_list = []
    for record in SeqIO.parse(monomers_path, "fasta"):
        monomers_list.append(record)
    return monomers_list

def get_cnt(cntID, str_lst_cen, str_lst_cnt):
    lst_cen = str_lst_cen.strip('][').split(', ')
    lst_cnt = str_lst_cnt.strip('][').split(', ')
    cnt = 0
    for i in range(len(lst_cen)):
        if cntID in lst_cen[i]:
            cnt += int(lst_cnt[i])
    return cnt

def main():
    args = parse_args()
    mnpath = os.path.join(args.indir, "hyb", "hybrid_mn.fasta")
    sumpath = os.path.join(args.indir, "summary_upgrade.csv")
    monomers_list = init_monomers(mnpath)

    cenLs = ["cen" + str(i) + "_" for i in range(1, 23)] + ["cenX_"]
    for cenID in cenLs:
        outfl = os.path.join(args.outdir, cenID + "mn.fa")
        outsm = os.path.join(args.outdir, cenID + "sum_up.csv")
        osf = open(outsm, "w")
        writer = csv.writer(osf)
        id = 0
        sum_upgr_full = []
        with open(outfl, "w") as fw:
            with open(sumpath) as sf:
                reader = csv.reader(sf)
                for line in reader:
                    if line[0] == "IterationId":
                        writer.writerow(["ID", "Name", "Final#"] + line)

                    if cenID in line[-6]:
                        sum_upgr_full.append([str(id), monomers_list[int(line[0])].id, get_cnt(cenID, line[-6], line[-5])] + line)
                        SeqIO.write(monomers_list[int(line[0])], fw, "fasta")
                        id += 1
            sum_upgr_full.sort(key=lambda x: -x[2])
            i = 0
            for line in sum_upgr_full:
                writer.writerow([i] + line)
                i += 1

        osf.close()


if __name__ == "__main__":
    main()