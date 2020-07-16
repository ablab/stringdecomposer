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

import edlib

def cnt_edist(lst):
    if len(str(lst[0])) == 0:
        return -1
    if len(str(lst[1])) == 0:
        return -1
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="locations")
    if result["editDistance"] == -1:
        return -1
    return 100 - result["editDistance"]*100//max(len(lst[0]), len(lst[1]))

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

def add_rc_monomers(monomers):
    res = []
    for m in monomers:
        res.append(m)
        res.append(make_record(m.seq.reverse_complement(), m.name + "'", m.id + "'"))
    return res

def convert_read(decomposition, read, monomers):
    res = []
    for d in decomposition:
        monomer, start, end = d["m"], d["start"], d["end"]
        scores = {}
        for m in monomers:
            score = cnt_edist([read.seq[start:end + 1], m.seq])
            scores[m.name] = score
        
        if monomer == None:
            for s in scores:
                if monomer == None or scores[s] > scores[monomer]:
                    monomer = s
        secondbest, secondbest_score = None, -1 
        for m in scores:
            if m != monomer and abs(scores[m] - scores[monomer]) < 5:
                if not secondbest or secondbest_score < scores[m]:
                    secondbest, secondbest_score = m, scores[m]
        res.append({"m": monomer, "start": str(d["start"]), "end": str(d["end"]), "score": scores[monomer], \
                                "second_best": str(secondbest), "second_best_score": secondbest_score, "alt": scores, "q": "+"})

    window = 2
    for i in range(len(res)):
        sm, cnt = 0, 0
        for j in range(i - window, i + window + 1):
            if j >= 0 and j < len(res):
                sm += res[j]["score"]
                cnt += 1
        if sm/cnt < 80:
            res[i]["q"] = "?"

    return res

def print_read(fout, fout_alt, dec, read, monomers):
    dec = convert_read(dec, read, monomers)
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

    reads = load_fasta(args.sequences, "map")
    monomers = load_fasta(args.monomers)
    monomers = add_rc_monomers(monomers)
    if args.decomposition.endswith("tsv"):
        convert_tsv(args.decomposition, reads, monomers, outfile)
    else:
        convert_fasta(args.decomposition, reads, monomers, outfile)

    
    

