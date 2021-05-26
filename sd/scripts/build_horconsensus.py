from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline

from Bio import AlignIO

import sys
import os
from os import listdir
from os.path import isfile, isdir, join
import argparse
import pandas as pd
import numpy as np

import edlib

MONOIDNT = 95
PATH = "/Sid/tdvorkina/monomers/hprc_sdruns/chm13-asm/"

hor_order = {"1": "T1", "2": "A2", "3": "R3", "4": "A4", "5": "A1/5/16/19", "6": "A6",\
             "7": "A7", "8": "A8", "9": "A4/9", "10": "A10", "11": "A11", "12": "H12",\
             "13": "K13/21", "14": "H14", "15": "K15", "16": "A16", "17":["R17", "A17+AO17"],\
             "18": "A18/20", \
             "19": "B19+D1/5/16/19+C19+G19", "20": "A20", "21": "K13/21", "22": "H14/22",\
              "X": "LX"}
hor_rev = {"3", "12", "13", "14", "15", "21", "22", "X"}

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

def load_monodec(cen):
    dec = []
    monomers = set()
    filename = os.path.join(PATH, "cen" + cen + "_dec.tsv")
    with open(filename, "r") as fin:
         for ln in fin.readlines():
             if len(ln.strip().split("\t")) < 5:
                 continue
             ref, mon, start, end, idnt = ln.strip().split("\t")[:5]
             dec.append([ref, mon, start, end, idnt])
             monomers.add(mon)
    #dec = revert_rcreads(dec)
    #dec = remove_badreads(dec)
    return dec, monomers


def shift(hor_lst, cen):
    min_ind = 0
    for i in range(len(hor_lst)):
        #if hor_lst[min_ind] > hor_lst[i]:
        if hor_lst[i] in hor_order[cen]:
            min_ind = i
            break
    return hor_lst[min_ind:] + hor_lst[:min_ind]

def load_horascycle(filename, monomers, cen):
    hors = []
    hor_name= ""
    mono_mp = {}
    cnt, hor_cnt = 0, 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if len(ln.split("\t")) < 2:
                continue
            hor_name, hor_seq = ln.strip().split("\t")[:-1]
            hor_lst = hor_seq.split(",")
            hor_lst = shift(hor_lst, cen)
            print(hor_lst)
            hor_nseq = ""
            for m in hor_lst:
                hor_nseq += str(monomers[m].seq)
            hors.append([hor_lst, hor_nseq])
    return hors

def run_clustal(mappings, pair_name):
    if len(mappings) == 1:
        save_fasta(PATH + "/pair_alns/" + pair_name + "_seq.fasta", mappings + mappings)
    else:
        save_fasta(PATH + "/pair_alns/" + pair_name + "_seq.fasta", mappings)

    clustalw_exe = r"/home/tdvorkina/soft/clustalo-1.2.4-Ubuntu-x86_64"
    #cline = ClustalwCommandline(clustalw_exe, infile=PATH + "/pair_alns/" + pair_name + "_seq.fasta", outfile=PATH + "/pair_alns/" + pair_name + "_seq.clu")
    #stdout, stderr = cline()

    filename = PATH + "/pair_alns/" + pair_name + "_seq.clu"
    total_alns = []
    with open(filename, "r") as fin:
         alns = []
         for ln in fin.readlines():
             ln = ln.replace("  ", " ")
             seq = "*"
             if len(ln.split()) >= 2:
                 seq = ln.strip().split()[1]
                 if seq[0] in {"A","C","G", "T", "-"}:
                     alns.append(seq)
             elif len(alns) > 0:
                  if len(total_alns) == 0 or len(total_alns) == len(alns):
                      if len(total_alns) == 0:
                          total_alns = ["" for _ in range(len(alns))]
                      for i in range(len(alns)):
                          total_alns[i] += alns[i]
                      alns = []
                  else:
                      print("Something went wrong")
                      exit(-1)
    return total_alns

def extract_consensus(total_alns):
    consensus = ""
    #print(total_alns[0])
    #print(total_alns[1])
    #print(total_alns[2])
    for i in range(len(total_alns[0])):
        score = {"A": 0, "C": 0, "G": 0, "T": 0, "-": 0}
        for j in range(len(total_alns)):
            score[total_alns[j][i]] += 1
        scores_lst = sorted([[it, score[it]] for it in score ], key = lambda x: -x[1])
        max_score = scores_lst[0][0]
        if scores_lst[0][1] == scores_lst[1][1]:
            max_score = "N"
        if max_score != "-":
            consensus += max_score
    return consensus

def align_mappings(mappings, pair_name):
    pair_name = pair_name.replace("/", "_")
    pair_name = pair_name.replace("(", "_")
    pair_name = pair_name.replace(")", "_")
    algns = run_clustal(mappings, pair_name)
    consensus = extract_consensus(algns)
    print(len(consensus))
    return consensus

def edist(lst):
    if len(str(lst[0])) == 0:
        return 100500
    if len(str(lst[1])) == 0:
        return 100500
    result = edlib.align(str(lst[0]), str(lst[1]), mode="SHW", task="path", k=100)
    if result["editDistance"] == -1:
        return 100500, []
    aln = edlib.getNiceAlignment(result, str(lst[0]), str(lst[1]))
    return result["editDistance"], aln

def glue_pairs(p1, p2):
    max_len = 200
    #print(p1)
    #print(" ".join(["" for _ in range(170)]), p2)
    eds = []
    ed, aln = edist([p1[-max_len:], p2])
    longest, longest_ind = 0, -1
    cur_len = 0
    #print(ed)
    #print(aln["query_aligned"])
    #print(aln["matched_aligned"])
    #print(aln["target_aligned"])
    for i in range(len(aln["matched_aligned"])):
        if aln["matched_aligned"][i] == "|":
            cur_len += 1
        else:
            if cur_len > longest:
                longest, longest_ind = cur_len, i - cur_len
            cur_len = 0
    if cur_len > longest:
        longest, longest_ind = cur_len, len(aln["matched_aligned"]) - cur_len
    #print(longest)
    #print(aln["query_aligned"][longest_ind: longest_ind + longest])
    #print(aln["target_aligned"][longest_ind: longest_ind + longest])
    i, j = len(p1) - max_len + len(aln["query_aligned"][:longest_ind].replace("-", "")), len(aln["target_aligned"][:longest_ind].replace("-", ""))
    #for ln in range(max_len, 100, -1):
    #    ed, aln = edist([p1[-ln:], p2])
    #    if ed < 50 and ln >= 150:
    #        print(ln, ed)
    #        print(aln["query_aligned"])
    #        print(aln["matched_aligned"])
    #        print(aln["target_aligned"])
    #    for i in range(len(p1)-ln):
    #        for j in range(len(p2) - ln):
    #            if p1[i: i + ln] == p2[j: j + ln]:
    #               print(i, j, ed, ln, len(p1[:i] + p2[j:]))
    #               return i, j
    print(i, j, ed, len(p1[:i] + p2[j:]))
    print("")
    return i, j

def build_consensus(hors, mono_dec, ref, cen):
    hors_seq = []
    for horfull in hors:
        hor = horfull[0]
        pairs_consensus = []
        #mappings = [make_record(Seq(horfull[1]), "fullHOR", "fullHOR")]
        for i in range(len(hor)):
            m1, m2 = hor[i], hor[(i+1)%len(hor)]
            mappings = []
            for j in range(len(monodec)-1):
                if monodec[j][1] == m1 and monodec[j+1][1] == m2 and float(monodec[j][4]) > 95 and float(monodec[j+1][4]) > 95:
                    s1, e1, s2, e2 = int(monodec[j][2]), int(monodec[j][3]), int(monodec[j+1][2]), int(monodec[j+1][3])
                    if s2 - e1 < 10:
                        mappings.append(make_record(ref[monodec[j][0]].seq[s1: e2 + 1], m1+"_" + m2 + str(s1), m1 + "_" +m2 + str(s1) ))
            print("pair: ", m1, m2, len(mappings))
            pairs_consensus.append(align_mappings(mappings, m1+"_" + m2))
        #align_mappings(mappings, cen + "_full")
        #exit()
        cur_consensus = ""
        border = []
        for i in range(len(pairs_consensus)):
            print("Pair ", i)
            border.append(glue_pairs(pairs_consensus[i], pairs_consensus[(i+1)%len(pairs_consensus)]))
        l, r = border[len(pairs_consensus)-1][1], 0
        for i in range(len(pairs_consensus)):
            r = border[i][0]
            if l > r:
               print("Something went wrong!", i, l, r)
               exit(-1)
            cur_consensus += pairs_consensus[i][l: r]
            l = border[i][1]
        print(len(cur_consensus))
        hors_seq.append(make_record(Seq(cur_consensus), "D" + cen + "Z1", "D" + cen + "Z1", str(len(cur_consensus)) + "bp " + ",".join(hor) ))
    return hors_seq

if __name__ == "__main__":
    cen = sys.argv[1]
    ref = load_fasta(os.path.join("/Bmo/kolga/data/SD/allCt", "cen" + cen + "ct.fa"), "map")
    monomers_seq = load_fasta("/Bmo/kolga/data/SD/HORmon/ValuableMonomers/mon" + cen + ".fa", "map")
    monodec, monomers = load_monodec(cen)
    hors = load_horascycle("/Sid/tdvorkina/monomers/HORmon/HORs/HOR" + cen + ".tsv", monomers_seq, cen)
    outfilename = "/Sid/tdvorkina/monomers/HORmon/HORs/HOR" + cen + ".fasta"
    consensus = build_consensus(hors, monodec, ref, cen)
    if cen in hor_rev:
        consensus_rev = []
        for c in consensus:
            consensus_rev.append(make_record(c.seq.reverse_complement(), c.name, c.id, c.description.split()[0] + " " + ",".join([d + "'" for d  in c.description.split()[1].split(",")[::-1] ]) ) )
        consensus = consensus_rev
    save_fasta(outfilename, consensus)
