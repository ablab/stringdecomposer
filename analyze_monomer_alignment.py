#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


import edlib

ED_THRESHOLD = 0.5

def cnt_edist(lst):
    if len(str(lst[0])) == 0:
        return -1
    if len(str(lst[1])) == 0:
        return -1
    ed_er = int(ED_THRESHOLD*len(lst[0]))
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path", k = ed_er)
    if result["editDistance"] == -1:
        return -1, None
    niceAlign = edlib.getNiceAlignment(result, str(lst[0]), str(lst[1]))
    return 100 - result["editDistance"]*100//max(len(lst[0]), len(lst[1])), niceAlign["query_aligned"], niceAlign["matched_aligned"], niceAlign["target_aligned"]

def cnt_pairwise(lst):
    if len(str(lst[0])) == 0:
        return -1
    if len(str(lst[1])) == 0:
        return -1
    alignment = pairwise2.align.globalms(str(lst[0]), str(lst[1]), 2, -1, -1, -.5)
    for i in range(1):
        print(alignment[i][0])
        print(alignment[i][1])
        print(alignment[i][2])
    

def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    else:
        records = list(SeqIO.parse(filename, "fasta"))
    return records

def load_decomposition_tsv(filename):
    reads_mapping = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt, q = ln.strip().split("\t")[:6]
            if sseqid not in reads_mapping:
                    reads_mapping[sseqid] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            rev = False
            if qseqid.endswith("'"):
                rev = True
                qseqid = qseqid[:-1]
            #if q == "+":
            reads_mapping[sseqid].append({"qid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt})
    for r in reads_mapping:
        if len(reads_mapping[r]) > 0 and reads_mapping[r][0]["rev"]:
            reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (-x["e"], x["s"]))
        else:
            reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (x["s"], -x["e"]))

    return reads_mapping

def print_alignments(hits, reads, monomers, prefix):
    with open(prefix + "monomers_consensus.fasta", "w") as fout4:
        for m in monomers:
            print(monomers[m].name)
            with open(prefix + monomers[m].name + ".fasta", "w") as fout:
                with open(prefix + monomers[m].name + "_edlib.tsv", "w") as fout2:
                    with open(prefix + monomers[m].name + "_consensus.tsv", "w") as fout3:
                        consensus = [{"A": 0, "C": 0, "G": 0, "T": 0, "": {}, "-": 0} for i in range(len(monomers[m].seq) + 1)]
                        num = 0
                        for r in hits:
                            read = reads[r]
                            for h in hits[r]:
                                if h["qid"] == monomers[m].name:
                                    if h["rev"]:
                                        alignment = read.seq[h["s"]:h["e"] + 1].reverse_complement()
                                    else:
                                        alignment = read.seq[h["s"]:h["e"] + 1]
                                    fout.write(">" + monomers[m].name + "_" + str(h["s"]) + "_" + str(h["e"]) + "_" + r[:10]  + "\n")
                                    fout.write(str(alignment) + "\n")
                                    score, query, align, target = cnt_edist([monomers[m].seq, alignment])
                                    cnt_pairwise([monomers[m].seq, alignment])
                                    ind = 1
                                    gap_str =""
                                    for i in range(len(query)):
                                        if align[i] == "|" or align[i] == ".":
                                            gap_str = ""
                                            consensus[ind][target[i]] += 1
                                            ind += 1
                                        elif align[i] == "-" and query[i] == "-":
                                            if len(gap_str) > 0:
                                                consensus[ind - 1][""][gap_str] -= 1  
                                            gap_str += target[i]
                                            if gap_str not in consensus[ind-1][""]:
                                                consensus[ind - 1][""][gap_str] = 0 
                                            consensus[ind-1][""][gap_str] += 1
                                        else:
                                            gap_str = ""
                                            consensus[ind]["-"] += 1
                                            ind += 1
                                    fout2.write( query + "\n" + "\t".join([align, str(len(align)), str(h["s"]), str(h["e"]), r[:10] ])+ "\n")
                                    num += 1
                            exit(-1)
                        print(consensus[0])
                        consensus_str = ""
                        keys = ["A", 'C', "G", "T", "-"]
                        for i in range(len(consensus)):
                            best_k = "A"
                            for k in keys:
                                if consensus[i][k] > consensus[i][best_k]:
                                    best_k = k
                            if best_k != "-" and consensus[i][best_k] > 0.6*num:
                                consensus_str += best_k
                            elif best_k != "-" and consensus[i][best_k] > 0.5*num:
                                consensus_str += best_k.lower()
                            elif best_k != "-" and consensus[i][best_k] <= 0.5*num:
                                print(["WTF ", consensus[i]])
                            best_gap = ""
                            if len(consensus[i][""]) != 0:
                                for s in consensus[i][""]:
                                    if len(best_gap) == 0  or consensus[i][""][s] > consensus[i][""][best_gap]:
                                        best_gap = s
                                if consensus[i][""][best_gap] > 0.6*num:
                                    consensus_str += best_gap
                                elif consensus[i][""][best_gap] > 0.5*num:
                                    consensus_str += best_gap.lower()
                        print(consensus_str)

                        fout3.write("\t".join(list(" -" + str(monomers[m].seq))) +"\n")
                        
                        for k in keys:
                            ss = [k]
                            for i in range(len(consensus)):
                                ss.append(str(consensus[i][k]))
                            fout3.write("\t".join(ss) + "\n")
                        ss = ["gap"]
                        for i in range(len(consensus)):
                            best_gap = ""
                            if len(consensus[i][""]) == 0:
                                ss.append("-")
                            else:
                                for s in consensus[i][""]:
                                    if len(best_gap) == 0  or consensus[i][""][s] > consensus[i][""][best_gap]:
                                        best_gap = s
                                ss.append(best_gap + "(" + str(consensus[i][""][best_gap]) + ")")
                        fout3.write("\t".join(ss) + "\n")
                        score, query, align, target = cnt_edist([monomers[m].seq, consensus_str.upper()])
                        fout3.write(query + " \n")
                        fout3.write(align + "\n")
                        fout3.write(target + " -- Consensus\n")
                        fout3.write(consensus_str + "\n")


            fout4.write(">" + monomers[m].name + "\n")
            fout4.write(consensus_str.upper() + "\n")            


# reads = load_fasta("../chr6/chm13_d6z1.rds_ge100kb.NHGRI.guppy_3.1.5_fixed.fa", "map")
# monomers = load_fasta("../chr6/D6Z1_hg38_monomers.fasta", "map")
# prefix = "../chr6/monomers_alignments/"
# hits = load_decomposition_tsv("../chr6/decomposition_reads_new.tsv")

# reads = load_fasta("../chrX/cenX_v0_8_2/final_sequence_4.fasta", "map")
# monomers = load_fasta("../chrX/X02418/monomers.fasta", "map")
# hits = load_decomposition_tsv("../chrX/decomposition_ref_new.tsv")
# prefix = "../chrX/monomers_alignments/"
# monomers = load_fasta("../chrX/DXZ1_rc_star_monomers.fasta", "map")
# hits = load_decomposition_tsv("../chrX/decomposition_ref_star_new.tsv")
# prefix = "../chrX/monomers_star_alignments/"
# reads = load_fasta("../chrX/results/centromeric_reads/centromeric_reads.fasta", "map")
# monomers = load_fasta("../chrX/X02418/monomers.fasta", "map")
# hits = load_decomposition_tsv("../chrX/decomposition_reads_new.tsv")
# prefix = "../chrX/monomers_alignments_reads/"

reads = load_fasta("../chrX/results/centromeric_reads/centromeric_reads.fasta", "map")
monomers = load_fasta("../chrX/DXZ1_rc_star_monomers.fasta", "map")
hits = load_decomposition_tsv("../chrX/decomposition_reads_star_new.tsv")
prefix = "./test_dir"
print_alignments(hits, reads, monomers, prefix)

