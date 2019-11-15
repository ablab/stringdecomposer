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
LOWEST_IDENTITY = 70

def cnt_edist(lst):
    if len(str(lst[0])) == 0:
        return -1
    if len(str(lst[1])) == 0:
        return -1
    ed_er = int(ED_THRESHOLD*len(lst[0]))
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="locations", k = ed_er)
    if result["editDistance"] == -1:
        return -1
    return 100 - result["editDistance"]*100//len(lst[0])

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
        if identity > LOWEST_IDENTITY:
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
                res.append([name, ind, m[2], m[0], m[1], m[3], idnt, "+"])
    new_res = []
    for i in range(len(res)):
        add = True
        for j in range(max(0, i - 3), i):
            if (res[i][1] + res[i][4]) - (res[j][1] + res[j][4]) < 50 or (res[i][1] + res[i][3]) - (res[j][1] + res[j][3]) < 50:
                if res[j][6] >= res[i][6]:
                    add = False

        for j in range(i + 1, min(i + 4, len(res))):
            if (res[j][1] + res[j][4]) - (res[i][1] + res[i][4]) < 50 or (res[j][1] + res[j][3]) - (res[i][1] + res[i][3]) < 50:
                if res[j][6] > res[i][6]:
                    add = False
        if add:
            new_res.append(res[i])

    WINDOW = 5
    idnts = []
    l = 0
    for it in new_res:
        idnts.append(it[6])
        sm = sum(idnts[l:])/(len(idnts) - l)
        if sm < 80:
            it[7] = "?"
        if len(idnts) > WINDOW:
            l += 1
    return new_res


def get_monomer_sequence(alns, avg_len, encoding):
    res = []
    prev_start, prev_end = -1, -1
    rc_num, total_len = 0, 0
    reverse = False
    for a in alns:
        name = a[2]
        ind = a[1]
        start, end = ind + a[3], ind + a[4]
        num_n = 0
        if prev_end != -1 and start - prev_end > 50:
            num_n = (start - prev_end)//avg_len
            if (start - prev_end)%avg_len*2 >= avg_len:
                num_n += 1
        for _ in range(num_n):
            res.append("?")
            total_len += 1
        res.append(encoding[name])
        total_len += 1
        if encoding[name].endswith("'"):
            rc_num += 1
        prev_start, prev_end = start, end
    if rc_num*2 > total_len:
        reverse = True
        for i in range(len(res)):
            if res[i].endswith("'"):
                res[i] = res[i][:-1]
            elif res[i] != "?":
                res[i] = res[i] + "'"
            else:
                res[i] = res[i]
        res = res[::-1]
    return "".join(res), "" if not reverse else "'", alns[0][1] + alns[0][3], alns[-1][1] + alns[-1][4]

def encode_monomers(monomers):
    if len(monomers) <= 26:
        res = {}
        i = "A"
        for m in monomers:
            if not m.name.endswith("'"):
                res[m.name] = i
                i = chr(ord(i) + 1)
            else:
                res[m.name] = res[m.name[:-1]] + "'"
        return res
    else:
        print("Warning: No encoding! Fasta won't be generated")
        return {}

def avg_length(monomers):
    res = 0
    for m in monomers:
        res += len(m.seq)
    return res//len(monomers)

def parallel_edlib_version(reads, monomers, outfile, t, identity_dif):
    LEN_STEP = 5000
    THREADS = int(t)
    SAVE_STEP = 300
    save_step = []
    new_reads = []
    avg_len = avg_length(monomers)
    encoding = encode_monomers(monomers)
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
    with open(outfile[:-len(".tsv")] + ".fasta", "w") as fout:
        fout.write("")
    
    start = 0
    for j in range(0, len(save_step)):
        print("Read " + new_reads[start].description + " aligning")
        all_ans = Parallel(n_jobs=THREADS)(delayed(slow_edlib_version)([new_reads[i], monomers, identity_dif]) for i in range(start, min(start + save_step[j], len(new_reads)) ))
        all_ans = transform_alignments(all_ans, new_reads, start)
        with open(outfile[:-len(".tsv")] + "_alt.tsv", "a+") as fout_alt:
            with open(outfile, "a+") as fout:
                for a in all_ans:
                    name = a[0]
                    ind = a[1]
                    second_best = 1 if str(a[5][0][0]) == str(a[2]) else 0
                    for alt in a[5]:
                        if str(alt[0]) == str(a[2]):
                            fout_alt.write("\t".join([name, str(alt[0]), str(ind + a[3]), str(ind + a[4]), "{:.2f}".format(alt[1]), "*"]) + "\n")
                        else:
                            fout_alt.write("\t".join([name, str(alt[0]), str(ind + a[3]), str(ind + a[4]), "{:.2f}".format(alt[1])]) + "\n")
                    if a[6] - a[5][second_best][1] < 5:
                        second_monomer, second_identity = a[5][second_best]
                    else:
                        second_monomer, second_identity = "None", -1
                    fout.write("\t".join([name, str(a[2]), str(ind + a[3]), str(ind + a[4]), "{:.2f}".format(a[6]), \
                                                                             second_monomer, "{:.2f}".format(second_identity)]) + "\n")
        if len(encoding) > 0 and len(all_ans) > 0:
            with open(outfile[:-len(".tsv")] + ".fasta", "a+") as fout:
                s, rev, left, right = get_monomer_sequence(all_ans, avg_len, encoding)
                fout.write(">" + new_reads[start].description  + rev + "/" + str(left) + "_" + str(right) + "\n")
                fout.write(s + "\n")
        print("Read " + new_reads[start].description + " aligned")
        start += save_step[j]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find highest score decomposition of each read into monomers')
    parser.add_argument('-s', '--sequences', help='fasta-file with long reads sequences', required=True)
    parser.add_argument('-m', '--monomers', help='fasta-file with monomers', required=True)
    parser.add_argument('-o', '--out',  help='output tsv-file, by default will be saved into decomposition.tsv', required=False)
    parser.add_argument('-i', '--identity', help='difference in identity for printed alignments, default i = 10', required=False)
    parser.add_argument('-t', '--threads', help='threads number', required=False)
    parser.add_argument('-a', '--stats', help='prints coverage and identity statistics for computed alignments', required=False, action='store_true')
    parser.add_argument('-d', '--identitythreshold', help='identity threshold for reliable mononomer alignment, default 70', required=False)

    args = parser.parse_args()
    t = args.threads
    if t == None:
        t = "1"
    i = args.identity
    if i == None:
        i = 10
    else:
        i = int(i)
    i = 100
    if args.identitythreshold == None:
        LOWEST_IDENTITY = 70
    else:
        LOWEST_IDENTITY = int(args.identitythreshold)

    print("Number of threads: " + t)
    outfile = args.out
    if outfile == None:
        outfile = "./decomposition.tsv"

    reads = load_fasta(args.sequences)
    monomers = load_fasta(args.monomers)
    monomers = add_rc_monomers(monomers)

    #if not os.path.exists(outfile):
    parallel_edlib_version(reads, monomers, outfile, t, i)

    if args.stats:
       sc = StatisticsCounter(outfile, args.sequences, args.monomers, LOWEST_IDENTITY)
       sc.save_stats(outfile[:-len(".tsv")])

