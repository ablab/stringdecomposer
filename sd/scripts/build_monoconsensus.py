from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqIO import FastaIO
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

import subprocess

import edlib
import random

MONOIDNT = 95
INF=10000000
clustalw_exe = r"/home/tdvorkina/soft/clustalo-1.2.4-Ubuntu-x86_64"

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
        #SeqIO.write(orfs, output_handle, "fasta")
        fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_out.write_file(orfs)

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

def load_bedfile(bedfile):
    dec = []
    monomers = set()
    with open(bedfile, "r") as fin:
         for ln in fin.readlines():
             if len(ln.strip().split("\t")) < 6:
                 continue
             ref, start, end, mon, idnt, rev = ln.strip().split("\t")[:6]
             dec.append([ref, mon, start, end, idnt, rev])
             monomers.add(mon)
    return dec, monomers

def shift(hor_lst):
    min_ind = 0
    for i in range(len(hor_lst)):
        if hor_lst[min_ind] > hor_lst[i]:
            min_ind = i
            break
    return hor_lst[min_ind:] + hor_lst[:min_ind]

def load_horascycle(filename):
    hors = []
    hor_name= ""
    mono_mp = {}
    cnt, hor_cnt = 0, 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if len(ln.split("\t")) < 2:
                continue
            hor_name, hor_seq = ln.strip().split("\t")[:2]
            hor_lst = hor_seq.split(",")
            #hor_lst = shift(hor_lst)
            print(hor_lst)
            hors.append([hor_name, hor_lst])
    return hors

def run_clustal(mappings, clustal_dir, pair_name):
    if not os.path.isdir(clustal_dir):
       os.makedirs(clustal_dir)
    if len(mappings) == 1:
        save_fasta(os.path.join(clustal_dir, pair_name + "_seq.fasta"), mappings + mappings)
    else:
        save_fasta(os.path.join(clustal_dir, pair_name + "_seq.fasta"), mappings)

    cline = ClustalwCommandline(clustalw_exe, infile=os.path.join(clustal_dir, pair_name + "_seq.fasta"), outfile=os.path.join(clustal_dir, pair_name + "_seq.clu"))
    stdout, stderr = cline()

    filename = os.path.join(clustal_dir, pair_name + "_seq.clu")
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
         if len(alns) > 0:
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
    scores = []
    #print(len(total_alns[0]), total_alns[0])
    #print(len(total_alns[1]), total_alns[1])
    #print(len(total_alns[2]), total_alns[2])
    for i in range(len(total_alns[0])):
        score = {"A": 0, "C": 0, "G": 0, "T": 0, "-": 0}
        for j in range(len(total_alns)):
            score[total_alns[j][i]] += 1
        scores_lst = sorted([[it, score[it]] for it in score ], key = lambda x: -x[1])
        scores.append(scores_lst)
        max_score = scores_lst[0][0]
        if scores_lst[0][1] == scores_lst[1][1]:
            max_score = "N"
        if max_score != "-":
            consensus += max_score
    return consensus, scores

def align_mappings(mappings, clustal_dir, pair_name):
    pair_name = pair_name.replace("/", "_")
    pair_name = pair_name.replace("(", "_")
    pair_name = pair_name.replace(")", "_")
    pair_name = pair_name.replace(".", "_")
    algns = run_clustal(mappings, clustal_dir, pair_name)
    consensus, scores = extract_consensus(algns)
    print(len(consensus))
    return consensus, scores

def edist(lst):
    if len(str(lst[0])) == 0:
        return INF, []
    if len(str(lst[1])) == 0:
        return INF, []
    result = edlib.align(str(lst[0]), str(lst[1]), mode="SHW", task="path", k=100)
    if result["editDistance"] == -1:
        return INF, []
    aln = edlib.getNiceAlignment(result, str(lst[0]), str(lst[1]))
    return result["editDistance"], aln

def glue_pairs(p1, p2):
    max_len = 200
    eds = []
    ed, aln = edist([p1[-max_len:], p2])
    longest, longest_ind = 0, -1
    cur_len = 0
    for i in range(len(aln["matched_aligned"])):
        if aln["matched_aligned"][i] == "|":
            cur_len += 1
        else:
            if cur_len > longest:
                longest, longest_ind = cur_len, i - cur_len
            cur_len = 0
    if cur_len > longest:
        longest, longest_ind = cur_len, len(aln["matched_aligned"]) - cur_len
    i, j = len(p1) - max_len + len(aln["query_aligned"][:longest_ind].replace("-", "")), len(aln["target_aligned"][:longest_ind].replace("-", ""))
    print(i, j, ed, len(p1[:i] + p2[j:]))
    print("")
    return i, j

def build_monoconsensus(monodec, ref, step, clustal_dir):
    mappings = {}
    for i in range(len(monodec)):
        m_name, m_s, m_e, m_idnt = monodec[i][1], int(monodec[i][2]), int(monodec[i][3]), float(monodec[i][4])
        if m_idnt > MONOIDNT:
            if m_name not in mappings:
                mappings[m_name] = []
            if monodec[i][5] == "+":
                mappings[m_name].append(make_record(ref[monodec[i][0]].seq[m_s - step: m_e + 1 + step], m_name + str(m_s), m_name + str(m_s)))
            else:
                mappings[m_name].append(make_record(ref[monodec[i][0]].seq[m_s - step: m_e + 1 + step].reverse_complement(), m_name + str(m_s) + "_rev", m_name + str(m_s) + "_rev"))
    m_consensus = []
    m_scores = {}
    for m in mappings:
        print(m)
        consensus, scores = align_mappings(mappings[m], clustal_dir, m)
        consensus, scores = consensus[step : len(consensus)-step], scores[step: len(scores)-step]
        cur_consensus = make_record(Seq(consensus), m, m)
        m_consensus.append(cur_consensus)
        m_scores[m] = scores
    return m_consensus, m_scores

def build_pairconsensus(hors, monodec, ref, clustal_dir):
    hors_seq = []
    for hor_desc in hors:
        hor = hor_desc[1]
        pairs_consensus = []
        for i in range(len(hor)):
            m1, m2 = hor[i], hor[(i+1)%len(hor)]
            mappings = []
            for j in range(len(monodec)-1):
                if ((monodec[j][1] == m1 and monodec[j+1][1] == m2 and monodec[j][5] == monodec[j+1][5] == "+") or\
                    (monodec[j][1] == m2 and monodec[j+1][1] == m1 and monodec[j][5] == monodec[j+1][5] == "-"))  \
                    and float(monodec[j][4]) > MONOIDNT and float(monodec[j+1][4]) > MONOIDNT:
                    s1, e1, s2, e2 = int(monodec[j][2]), int(monodec[j][3]), int(monodec[j+1][2]), int(monodec[j+1][3])
                    if s2 - e1 < 10:
                        if monodec[j][5] == "+":
                            mappings.append(make_record(ref[monodec[j][0]].seq[s1: e2 + 1], m1+"_" + m2 + str(s1), m1 + "_" +m2 + str(s1) ))
                        else:
                            mappings.append(make_record(ref[monodec[j][0]].seq[s1: e2 + 1].reverse_complement(), m1+"_" + m2 + str(s1) + "_rev", m1 + "_" +m2 + str(s1) + "_rev" ))
            print("pair: ", m1, m2, len(mappings))
            if len(mappings) > 0:
                pairs_consensus.append(align_mappings(mappings, clustal_dir, m1+"_" + m2)[0])
        if len(pairs_consensus) > 0:
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
            hors_seq.append(make_record(Seq(cur_consensus), hor_desc[0], hor_desc[0], str(len(cur_consensus)) + "bp " + ",".join(hor) ))
    return hors_seq


def run(sequences, monomers, num_threads, scoring, batch_size, raw_file):
    ins, dels, mm, match = scoring.split(",")
    p = os.path.abspath(__file__)
    sd_exec_file = p[:-len("/scripts/build_monoconsensus.py")] + "/bin/dp"
    print("Run: ", sd_exec_file, sequences, monomers, num_threads, batch_size, scoring, file=sys.stderr)
    with open(raw_file, 'w') as f:
        subprocess.run([sd_exec_file, sequences, monomers, num_threads, batch_size, ins, dels, mm, match], stdout = f, check = True)
    with open(raw_file, 'r') as f:
        raw_decomposition = "".join(f.readlines())
    return raw_decomposition


def divide_into_monomers(hors_lst, hors, monomers, horfile, monofile, outtsv):
    hors_res, res, bed = [], [], []
    for i in range(len(hors)):
        h, start_mono, hor_len = hors[i], hors_lst[i][1][0], len(hors_lst[i][1])
        triple_hors = []
        triple_hors.append(make_record(h.seq + h.seq + h.seq, h.id, h.name, h.description))
        save_fasta(horfile, triple_hors)
        decomposition = run(horfile, monofile, "1", "-1,-1,-1,1", "5000", outtsv)
        inHOR, start_shift = False, 0
        newhor_seq = ""
        for ln in decomposition.split("\n")[1:-1]:
            ref, mono, start, end, score = ln.split("\t")[:5]
            if mono == start_mono:
                inHOR, start_shift = True, int(start)
                print("shift", start_shift)
            if inHOR and len(res) < hor_len:
                res.append(make_record(triple_hors[0].seq[int(start): int(end) + 1], mono, mono, ref + ":" + str(int(start) - start_shift) + "-" + str(int(end) + 1 - start_shift) ))
                r, g, b = random.randint(0, 255), random.randint(0, 255), random.randint(0, 255)
                bed.append("\t".join([ref, str(int(start) - start_shift), str(int(end) + 1 - start_shift), mono, str(int(float(score))), "+", str(int(start) - start_shift), str(int(end) + 1 - start_shift), ",".join([str(r), str(g), str(b)]) ]))
                newhor_seq += triple_hors[0].seq[int(start): int(end) + 1]
        hors_res.append(make_record(newhor_seq, h.id, h.id, str(len(newhor_seq)) + "bp " + ",".join(hors_lst[i][1]) ))
    return hors_res, res, bed


def main():
    parser = argparse.ArgumentParser(description='Extracts monomer/HOR consensus from annotation')
    parser.add_argument('sequences', help='fasta-file with annotated sequences')
    parser.add_argument('annotation', help='bed-file with annotation')
    parser.add_argument('outdir', help='output directory')
    parser.add_argument('--hors',  help='tsv-file with HOR description', required=False)
    parser.add_argument('--extend',  help='number of bp to extend monomer alignment', default=4, required=False)
    parser.add_argument('--seed',  help='seed for colors generation in bed-files', default=123, required=False)
    args = parser.parse_args()
    
    random.seed(int(args.seed))
    ref = load_fasta(args.sequences, "map")
    bedfile = args.annotation
    
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    mono_outfilename = os.path.join(args.outdir, "monomer_consensus.fasta")
    monodec, monomers = load_bedfile(bedfile)
    consensus, scores = build_monoconsensus(monodec, ref, int(args.extend), os.path.join(args.outdir, "clustal_alns"))
    save_fasta(mono_outfilename, consensus)
    with open(os.path.join(args.outdir, "monomer_consensus.tsv"), "w") as fout:
        for m in scores:
            for i in range(len(scores[m])):
                s = scores[m][i]
                fout.write(m + "\t" + str(i) + "\t" + "\t".join(x[0] + ":" + str(x[1]) for x in s) + "\n")
    print("Monomer consensus can be found", mono_outfilename)

    if args.hors != None:
       hors_tsv = args.hors
       hors = load_horascycle(hors_tsv)
       hor_consensus = build_pairconsensus(hors, monodec, ref, os.path.join(args.outdir, "clustal_alns"))
       if len(hor_consensus) > 0:
           hor_outfilename = os.path.join(args.outdir, "hor_consensus.fasta")
           hor_consensus_shifted, pair_monomers, bed = divide_into_monomers(hors, hor_consensus, consensus, hor_outfilename, mono_outfilename, os.path.join(args.outdir, "sd_raw.tsv"))
           hor_outfilename = os.path.join(args.outdir, "hor_consensus.fasta")
           save_fasta(hor_outfilename, hor_consensus_shifted)
           print("HOR consensus can be found", hor_outfilename)
           pmono_outfilename = os.path.join(args.outdir, "monomer_paired_consensus.fasta")
           save_fasta(pmono_outfilename, pair_monomers)
           print("Paired monomer consensus can be found", pmono_outfilename)
           outfilename = os.path.join(args.outdir, "monomer_paired_consensus.bed")
           with open(outfilename, "w") as fout:
               for ln in bed:
                   fout.write(ln + "\n")
           print("Paired monomer consensus bed can be found", outfilename)
       else:
           print("Paired wasn't generated - no pairs found")

if __name__ == "__main__":
    main()