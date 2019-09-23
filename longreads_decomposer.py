#!/usr/bin/env python
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import cProfile

import os
import sys
import argparse

import re
import numpy as np

from joblib import Parallel, delayed
import edlib

from monomers_alignment_statistics import StatisticsCounter


ED_THRESHOLD = 0.5

def cnt_edist(lst):
    if len(str(lst[0])) == 0:
        return -1
    if len(str(lst[1])) == 0:
        return -1
    ed_er = int(ED_THRESHOLD*len(lst[0]))
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="locations", k = ed_er)
    if result["editDistance"] == -1:
        return -1
    return 100 - result["editDistance"]*100//max(len(lst[0]), len(lst[1]))

def cnt_suffix_edist(lst):
    if len(str(lst[0])) == 0:
        return -1, -1
    if len(str(lst[1])) == 0:
        return -1, -1
    ed_er = int(ED_THRESHOLD*len(lst[0]))
    result = edlib.align(str(lst[0])[::-1], str(lst[1])[::-1], mode="SHW", task="locations", k = ed_er)
    if result["editDistance"] == -1:
        return -1, -1 
    return 100 - result["editDistance"]*100//len(lst[0]), len(lst[1]) - result["locations"][0][1] - 1

def load_fasta(filename):
    records = [rec.upper() for rec in SeqIO.parse(filename, "fasta")]
    return records

def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)

def add_rc_monomers(monomers):
    res = []
    for m in monomers:
        res.append(m)
        res.append(make_record(m.seq.reverse_complement(), m.name + "'", m.id + "'"))
    return res

def choose_best(read, read_dp, i, monomers, identity_dif):
    best_score, best_ind, best_monomer = -1, -1, ""
    identities = []
    for m in monomers:
        identity, start = cnt_suffix_edist([m.seq, read[: i + 1]])
        prev = 0
        if start - 1 >= 0:
            prev = read_dp[start-1]
        if identity > 70:
            if identity + prev > best_score:
                best_score, best_ind, best_monomer = identity + prev, start, m.name
    if best_score > -1:
        best_idnt = best_score 
        if best_score - read_dp[best_ind - 1] >= 0:
            best_idnt -= read_dp[best_ind - 1]
        for m in monomers:
            idnt = cnt_edist([m.seq, read[best_ind: i + 1] ])
            if idnt + identity_dif >= best_idnt and best_idnt != -1:
                identities.append([m.name, idnt])
        identities = sorted(identities, key = lambda x: -x[1])
    return best_score, best_ind, best_monomer, identities

def slow_edlib_version(args):
    r, monomers, identity_dif = args[0], args[1], args[2]
    read_len = len(r.seq)
    read_dp = [0 for _ in range(read_len)]
    ans_dp = [0 for _ in range(read_len)]
    monomer_dp = ["" for _ in range(read_len)]
    identity_dp = [[] for _ in range(read_len)]
    ans_dp[0] = -1
    for i in range(1, read_len):
        read_dp[i], ans_dp[i], monomer_dp[i], identity_dp[i] = choose_best(r.seq, read_dp, i, monomers, identity_dif)
        if read_dp[i] < read_dp[i - 1]:
            read_dp[i], ans_dp[i], monomer_dp[i], identity_dp[i] = read_dp[i-1], i - 1, "", []
    ans = []
    ind = read_len - 1
    while ind >= 0:
        if monomer_dp[ind] == "" and len(ans) > 0 and ans[-1][2] == "":
            ans[-1] = ans_dp[ind], ans[-1][1], "", ""
        else:
            ans.append([ans_dp[ind], ind, monomer_dp[ind], identity_dp[ind]])
        ind = ans_dp[ind]

    return [r.name, ans[::-1]]

def transform_alignments(alns, new_reads, s):
    res = []
    prev = {"name": None, "start": None, "end": None, "idnt": None}
    for i in range(len(alns)):
        a = alns[i]
        name = new_reads[s + i].description
        for m in a[1]:
            if m[0] != -1 and m[2] != "":
                ind = int(a[0])
                cur_name = m[2]
                start, end = ind + m[0], ind + m[1]
                for aa in m[3]:
                    if aa[0] == cur_name:
                        idnt = aa[1]
                        break
                if cur_name == prev["name"] and end - prev["end"] < 150:
                    if idnt > prev["idnt"]:
                        res[-1] = [name, ind, m[2], m[0], m[1], m[3], idnt, "+"]
                else:
                    res.append([name, ind, m[2], m[0], m[1], m[3], idnt, "+"])
                prev = {"name": cur_name, "start": start, "end": end, "idnt": idnt}

    WINDOW = 5
    idnts = []
    l = 0
    for it in res:
        idnts.append(it[6])
        sm = sum(idnts[l:])/(len(idnts) - l)
        if sm < 80:
            it[7] = "?"
        if len(idnts) > WINDOW:
            l += 1
    return res        


def parallel_edlib_version(reads, monomers, outfile, t, identity_dif):
    LEN_STEP = 5000
    THREADS = int(t)
    SAVE_STEP = 300
    save_step = []
    new_reads = []
    for r in reads:
        cnt = 0
        for i in range(0, len(r.seq), LEN_STEP):
            new_reads.append(make_record(r.seq[i: min(i + LEN_STEP + 200, len(r.seq))], str(i), str(i), r.name))
            cnt += 1
        save_step.append(cnt)
    print("Initial number of reads: " + str(len(reads)) + ", Divided into chunks and reverse complement: " + str(len(new_reads)))
    with open(outfile, "w") as fout:
        fout.write("")
    with open(outfile[:-len(".tsv")] + "_alt.tsv", "w") as fout:
        fout.write("")
    
    start = 0
    for j in range(0, len(save_step)):
        all_ans = Parallel(n_jobs=THREADS)(delayed(slow_edlib_version)([new_reads[i], monomers, identity_dif]) for i in range(start, min(start + save_step[j], len(new_reads)) ))
        all_ans = transform_alignments(all_ans, new_reads, start)
        print("Read " + new_reads[start].description + " aligned")
        with open(outfile[:-len(".tsv")] + "_alt.tsv", "a+") as fout_alt:
            with open(outfile, "a+") as fout:
                for a in all_ans:
                    name = a[0]
                    ind = a[1]
                    fout.write("\t".join([name, str(a[2]), str(ind + a[3]), str(ind + a[4]), "{:.2f}".format(a[6]), a[7]]) + "\n")
                    add_star = True if len(a[5]) > 1 else False 
                    for alt in a[5]:
                        if add_star:
                            fout_alt.write("\t".join([name, str(alt[0]), str(ind + a[3]), str(ind + a[4]), "{:.2f}".format(alt[1]), "*"]) + "\n")
                        else:
                            fout_alt.write("\t".join([name, str(alt[0]), str(ind + a[3]), str(ind + a[4]), "{:.2f}".format(alt[1])]) + "\n")
        start += save_step[j]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find highest score decomposition of each read into monomers')
    parser.add_argument('-s', '--sequences', help='fasta-file with long reads sequences', required=True)
    parser.add_argument('-m', '--monomers', help='fasta-file with monomers', required=True)
    parser.add_argument('-o', '--out',  help='output tsv-file, by default will be saved into decomposition.tsv', required=False)
    parser.add_argument('-i', '--identity', help='difference in identity for printed alignments, default i = 10', required=False)
    parser.add_argument('-t', '--threads', help='threads number', required=False)

    args = parser.parse_args()
    t = args.threads
    if t == None:
        t = "1"
    i = args.identity
    if i == None:
        i = 10
    else:
        i = int(i)
    print("Number of threads: " + t)
    outfile = args.out
    if outfile == None:
        outfile = "./decomposition.tsv"

    reads = load_fasta(args.sequences)
    monomers = load_fasta(args.monomers)

    monomers = add_rc_monomers(monomers)

    #cProfile.run("parallel_edlib_version(reads, monomers, outfile, t, i)")
    stats_only = True
    if not stats_only:
        parallel_edlib_version(reads, monomers, outfile, t, i)

    if stats_only:
       sc = StatisticsCounter(outfile, args.sequences, args.monomers)
       sc.save_stats(outfile[:-len(".tsv")])

