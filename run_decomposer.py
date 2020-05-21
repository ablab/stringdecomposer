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

import subprocess

import numpy as np; np.random.seed(0)
import pandas as pd

import re
import edlib

import joblib
from itertools import groupby

from nonmono_regions_extraction import contruct_nm_regions

p = os.path.abspath(__file__)
logreg_file = p[:-len("run_decomposer.py")] + "/models/new_ont_logreg_model.sav"
clf = joblib.load(logreg_file)

def edist(lst):
    if len(str(lst[0])) == 0:
        return -1, ""
    if len(str(lst[1])) == 0:
        return -1, ""
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path")
    return result["editDistance"], result["cigar"]

def aai(ar):
    p1, p2 = str(ar[0]), str(ar[1])
    if p1.endswith("*"):
        p1 = p1[:-1]
    if p2.endswith("*"):
        p2 = p2[:-1]
    ed, cigar = edist([str(p1), str(p2)])
    if ed == -1:
        return 0
    total_length = 0 #max(len(p1), len(p2))
    n = 0
    for c in cigar:
        if c.isdigit():
            n = n*10 + int(c)
        else:
            total_length += n
            n = 0
    matches = re.findall(r'\d+=', cigar)
    aai = 0.0
    for m in matches:
        aai += int(m[:-1])
    aai /= total_length
    return aai*100


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

def add_rc_monomers(monomers):
    res = []
    for m in monomers:
        res.append(m)
        res.append(make_record(m.seq.reverse_complement(), m.name + "'", m.id + "'"))
    return res

def convert_to_homo(seq):
    res = ''.join(i for i, _ in groupby(seq))
    return res

def classify(reads_mapping):
    df = pd.DataFrame(reads_mapping)
    df["idnt_diff"] = df["score"] - df["second_best_score"]
    X = pd.concat([df["score"], df["idnt_diff"]], axis=1, keys = ["idnt", "idnt_diff"])
    X_scaled = X
    y_pred = list(clf.predict(X_scaled))
    for i in range(len(reads_mapping)):
        if y_pred[i] != 1:
            reads_mapping[i]["q"] = "?" 
    return reads_mapping

def convert_read(decomposition, read, monomers):
    res = []
    for d in decomposition:
        monomer, start, end = d["m"], d["start"], d["end"]
        scores = {}
        for m in monomers:
            score = aai([read.seq[start:end + 1], m.seq])
            scores[m.name] = score
        
        if monomer == None:
            for s in scores:
                if monomer == None or scores[s] > scores[monomer]:
                    monomer = s
        secondbest, secondbest_score = None, -1
        for m in scores:
            if m != monomer: # and abs(scores[m] - scores[monomer]) < 5:
                if not secondbest or secondbest_score < scores[m]:
                    secondbest, secondbest_score = m, scores[m]

        homo_scores = []
        homo_subseq = convert_to_homo(read.seq[start:end + 1])
        for m in monomers:
            score = aai([homo_subseq, convert_to_homo(m.seq)])
            homo_scores.append([m.name, score])
        homo_scores = sorted(homo_scores, key = lambda x: -x[1])
        res.append({"m": monomer, "start": str(d["start"]), "end": str(d["end"]), "score": scores[monomer], \
                                "second_best": str(secondbest), "second_best_score": secondbest_score,\
                                "homo_best": homo_scores[0][0], "homo_best_score": homo_scores[0][1],\
                                "homo_second_best": homo_scores[1][0], "homo_second_best_score": homo_scores[1][1],\
                                "alt": scores, "q": "+"})

    res = classify(res)
    return res

def print_read(fout, fout_alt, dec, read, monomers, identity_th):
    dec = convert_read(dec, read, monomers)
    for d in dec:
        if d["score"] >= identity_th:
            fout.write("\t".join([read.name, d["m"], d["start"], d["end"], "{:.2f}".format(d["score"]), \
                                                    d["second_best"], "{:.2f}".format(d["second_best_score"]), \
                                                    d["homo_best"], "{:.2f}".format(d["homo_best_score"]), \
                                                    d["homo_second_best"], "{:.2f}".format(d["homo_second_best_score"]), d["q"]]) + "\n")
            for a in d["alt"]:
                star = "-"
                if a == d["m"]:
                    star = "*"
                fout_alt.write("\t".join([read.name, a, d["start"], d["end"], "{:.2f}".format(d["alt"][a]), star]) + "\n")

def convert_tsv(decomposition, reads, monomers, outfile, identity_th):
    with open(outfile[:-len(".tsv")] + "_alt.tsv", "w") as fout_alt:
        with open(outfile, "w") as fout:
            cur_dec = []
            prev_read = None
            for ln in decomposition.split("\n")[:-1]:
                read, monomer, start, end = ln.split("\t")[:4]
                read = read.split()[0]
                monomer = monomer.split()[0]
                if read != prev_read and prev_read != None:
                    print_read(fout, fout_alt, cur_dec, reads[prev_read], monomers, identity_th)
                    cur_dec = []
                prev_read = read
                start, end = int(start), int(end)
                cur_dec.append({"m": monomer, "start": start, "end": end})
            if len(cur_dec) > 0:
                print_read(fout, fout_alt, cur_dec, reads[prev_read], monomers, identity_th)

def form_nm_decomposition(non_mono, clusters, reads, outfile):
    non_mono_byreads = {}
    for m in non_mono:
        cl_id = int(m.id.split("_")[1])
        for region in clusters[cl_id]:
            if region["r"] not in non_mono_byreads:
                non_mono_byreads[region["r"]] = []
            if region["rev"]:
                idnt = aai([m.seq.reverse_complement(), reads[region["r"]].seq[region["s"]:region["e"] + 1] ])
            else:
                idnt = aai([m.seq, reads[region["r"]].seq[region["s"]:region["e"] + 1] ])
            non_mono_byreads[region["r"]].append({"m": m.id, "s": region["s"], "e": region["e"], "rev": region["rev"], "idnt": idnt})
    for r in non_mono_byreads:
        non_mono_byreads[r] = sorted(non_mono_byreads[r], key = lambda x: x["s"])
    decomposition = []
    with open(outfile[:-len(".tsv")] + "_with_masking.tsv", "w") as fout:
        with open(outfile, "r") as fin:
            cur_read, cur_ind = None, 0
            in_nonmono_region = False
            for ln in fin.readlines():
                read, monomer, start, end = ln.split("\t")[:4]
                read = read.split()[0]
                monomer = monomer.split()[0]
                if cur_read != read:
                    cur_read = read
                    cur_ind = 0
                start, end = int(start), int(end)
                if read in non_mono_byreads and cur_ind < len(non_mono_byreads[read]):
                    if start == non_mono_byreads[read][cur_ind]["s"]:
                        in_nonmono_region = True
                        name, n_start, n_end = non_mono_byreads[read][cur_ind]["m"], non_mono_byreads[read][cur_ind]["s"], non_mono_byreads[read][cur_ind]["e"]
                        if non_mono_byreads[read][cur_ind]["rev"]:
                            name += "'"
                        fout.write("\t".join([read, name, str(n_start), str(n_end), "{:.2f}".format(non_mono_byreads[read][cur_ind]["idnt"]), \
                                                    "None", "{:.2f}".format(-1), \
                                                    "None", "{:.2f}".format(-1), \
                                                    "None", "{:.2f}".format(-1), "+"]) + "\n")
                if not in_nonmono_region:
                    fout.write(ln)
                if read in non_mono_byreads and cur_ind < len(non_mono_byreads[read]):
                    if end == non_mono_byreads[read][cur_ind]["e"]:
                        in_nonmono_region = False
                        cur_ind += 1

def run(sequences, monomers, num_threads, scoring, batch_size, raw_file):
    ins, dels, mm, match = scoring.split(",")
    p = os.path.abspath(__file__)
    sd_exec_file = p[:-len("run_decomposer.py")] + "/src/dp"
    print("Run", sd_exec_file, " with parameters ", sequences, monomers, num_threads, batch_size, scoring, file=sys.stderr)
    with open(raw_file, 'w') as f:
        subprocess.run([sd_exec_file, sequences, monomers, num_threads, batch_size, ins, dels, mm, match], stdout = f, check = True)
    with open(raw_file, 'r') as f:
        raw_decomposition = "".join(f.readlines())
    return raw_decomposition

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Decomposes string into blocks alphabet')
    parser.add_argument('sequences', help='fasta-file with long reads or genomic sequences')
    parser.add_argument('monomers', help='fasta-file with monomers')
    parser.add_argument('-d', '--data-type',  help='type of reads (hifi or ont)', choices=["hifi", "ont"], required=True)
    parser.add_argument('-t', '--threads',  help='number of threads (by default 1)', default="1", required=False)
    parser.add_argument('-o', '--out-file',  help='output tsv-file (by default final_decomposition.tsv)', default="./final_decomposition.tsv", required=False)
    parser.add_argument('-i', '--min-identity',  \
                         help='only monomer alignments with percent identity >= MIN_IDENTITY are printed (by default MIN_IDENTITY=0)', type=int, default=0, required=False)
    parser.add_argument('-s', '--scoring', \
                         help='set scoring scheme for SD in the format "insertion,deletion,mismatch,match" (by default "-1,-1,-1,1")', default="-1,-1,-1,1", required=False)
    parser.add_argument('-b', '--batch-size',  help='set size of the batch in parallelization (by default 5000)', type=str, default="5000", required=False)
    parser.add_argument('-m', '--mask',  help='identifies similar large non-monomeric regions in reads and mask them as NM_X (where X is identifier of the region)', action="store_true")

    args = parser.parse_args()
    raw_decomposition = run(args.sequences, args.monomers, args.threads, args.scoring, args.batch_size, args.out_file[:-len(".tsv")] + "_raw.tsv")
    print("Saved raw decomposition to " + args.out_file[:-len(".tsv")] + "_raw.tsv", file=sys.stderr)

    reads = load_fasta(args.sequences, "map")
    monomers = load_fasta(args.monomers)
    monomers = add_rc_monomers(monomers)
    print("Transforming raw alignments...", file=sys.stderr)
    convert_tsv(raw_decomposition, reads, monomers, args.out_file, int(args.min_identity))
    print("Transformation finished. Results can be found in " + args.out_file, file=sys.stderr)

    if args.mask:
        print("Searching for non-monomeric regions started..")
        non_mono, clusters = contruct_nm_regions(args.out_file, reads, args.data_type)
        for m in non_mono:
            print(m.name)
        print(len(non_mono), len(clusters))
        new_monomer_file = "/".join(args.out_file.split("/")[:-1]) + "/monomers_with_nm_regions.fasta"
        new_dec_elements = [x for x in monomers if not x.id.endswith("'")] + non_mono
        print("Saving new set of elements to decompose to ", new_monomer_file)
        save_fasta(new_monomer_file, new_dec_elements)
        print("Saving new decomposition to ", args.out_file [:-len(".tsv")] + "_with_masking.tsv")
        form_nm_decomposition(non_mono, clusters, reads, args.out_file)
        #print("Rerunning SD on reads with non-monomeric regions and new set of elements to decompose..")
