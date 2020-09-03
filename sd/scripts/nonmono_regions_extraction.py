import sys
import edlib

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

import argparse
import os, shutil
import re

def edist(lst):
    if len(str(lst[0])) == 0:
        return 100500
    if len(str(lst[1])) == 0:
        return 100500
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW")
    return result["editDistance"]

def edist_hw(lst):
    if len(str(lst[0])) == 0:
        return 100500
    if len(str(lst[1])) == 0:
        return 100500
    result = edlib.align(str(lst[0]), str(lst[1]), mode="HW", task="path")
    niceAlign = edlib.getNiceAlignment(result, str(lst[0]), str(lst[1]))
    print(niceAlign["query_aligned"])
    print(niceAlign["matched_aligned"])
    print(niceAlign["target_aligned"])
    print(result["editDistance"])
    print("")
    return result["editDistance"]

def edist_aai(lst):
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
    ed, cigar = edist_aai([str(p1), str(p2)])
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


def load_regions(filename):
    reads_regions = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt  = ln.strip().split("\t")[:5]
            qual = ln.strip().split("\t")[-1]
            sseqid = sseqid.split()[0]
            rev = qseqid.endswith("'")
            if sseqid not in reads_regions:
                reads_regions[sseqid] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            reads_regions[sseqid].append({"idnt": idnt, "s": s, "e": e, "q": qual, "qid": qseqid, "rev": rev})
    for r in reads_regions:
        reads_regions[r] = sorted(reads_regions[r], key=lambda x: x["s"])
    return reads_regions

def extract_nonmono_regions(reads_regions, reads, datatype, th_border = 95, th_mono = 85):
    if datatype == "asm":
        window = 20
    else:
        window = 5
    for r in reads_regions:
        l, e = -1, -1
        idnt_sum = 0
        for j in range(len(reads_regions[r])):
            if j - window >= 0:
                idnt_sum -= reads_regions[r][j - window]["idnt"]
            idnt_sum += reads_regions[r][j]["idnt"]
            if l == -1 and idnt_sum/min(window, j + 1) >= th_border:
                l = j
            if idnt_sum/window >= th_border:
                e = j - window
        if e -l > 20:
            #print(reads_regions[r][l]["s"], reads_regions[r][e]["e"])
            reads_regions[r] = reads_regions[r][l:e + 1]
        else:
            reads_regions[r] = []

    collapsed_nonmono_regions = {}
    for r in reads_regions:
        if len(reads_regions[r]) > 0:
            in_region = False
            collapsed_nonmono_regions[r] = []
            cur_idnt = 0
            for i in range(len(reads_regions[r])):
                if reads_regions[r][i]["idnt"] < th_mono or reads_regions[r][i]["q"] == "?":
                    if in_region:
                        collapsed_nonmono_regions[r][-1].append(reads_regions[r][i])
                    else:
                        collapsed_nonmono_regions[r].append([reads_regions[r][i-1]])
                        in_region = True
                        #print("new region")
                        collapsed_nonmono_regions[r][-1].append(reads_regions[r][i])
                    #print(reads_regions[r][i]["s"], reads_regions[r][i]["e"],reads_regions[r][i]["idnt"], th_mono)
                else:
                    if in_region:
                        collapsed_nonmono_regions[r][-1].append(reads_regions[r][i])
                    in_region = False
    new_collapsed_nonmono_regions = {}
    for r in collapsed_nonmono_regions:
        if len(collapsed_nonmono_regions[r]) > 0:
            new_collapsed_nonmono_regions[r] = []
            for region in collapsed_nonmono_regions[r]:
                new_collapsed_nonmono_regions[r].append({ "seq": "".join(reads[r].seq[region[1]["s"] : region[-2]["e"] + 1]), \
                                                          "s": region[1]["s"], "e": region[-2]["e"], \
                                                          "prev_id": region[0]["qid"], "prev_seq": reads[r].seq[region[0]["s"]: region[0]["e"] + 1],
                                                          "next_id": region[-1]["qid"], "next_seq": reads[r].seq[region[-1]["s"]: region[-1]["e"] + 1]})
    return new_collapsed_nonmono_regions

def check_corrupted_as(seq, monomers):
    min_ed = len(seq)
    alsat = ""
    min_name = None
    for m in monomers:
        print(m.name)
        if len(m.seq) < len(seq):
            ed = edist_hw([m.seq, seq])
            min_ed = min(min_ed, min(ed, edist_hw([m.seq.reverse_complement(), seq])))
        else:
            ed = edist_hw([seq, m.seq])
            min_ed = min(min_ed, min(ed, edist_hw([seq, m.seq.reverse_complement()])))
    print(min_ed, len(seq))
    if min_ed < max(20, 0.3*len(seq)):
        alsat = "CM"
    else:
        print(len(seq), min_ed)
    return alsat

def restore_aln(dp, i, j, match, mismatch, gap, s1, s2):
    aln = [[],[]]
    x, y = i, j
    while x >= 0 and y >= 0 and dp[x][y] != 0:
        if (s1[y-1] == s2[x-1] and dp[x][y] == dp[x-1][y-1] + match) or (s1[y-1] != s2[x-1] and dp[x][y] == dp[x-1][y-1] + mismatch):
            aln[0].append(s1[y-1])
            aln[1].append(s2[x-1])
            x -= 1
            y -= 1
        elif dp[x][y] == dp[x-1][y] + gap:
            aln[0].append("-")
            aln[1].append(s2[x-1])
            x -= 1
        elif dp[x][y] == dp[x][y-1] + gap:
            aln[0].append(s1[y-1])
            aln[1].append("-")
            y -= 1
    return [aln[0][::-1], aln[1][::-1]]


def local_alignments(seq1, seq2, minlen = 20):
    seq1, seq2 = collapse_homo(seq1), collapse_homo(seq2)
    match, mismatch, gap = 2, -3, -2
    dp = [[0 for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    aln_len = [[0 for _ in range(len(seq1) + 1)] for _ in range(len(seq2) + 1)]
    alns = []
    for i in range(1, len(seq2) + 1):
        for j in range(1, len(seq1) + 1):
            if seq1[j - 1] == seq2[i - 1]:
                dp[i][j] = max(dp[i][j], dp[i-1][j-1] + match)
                if dp[i][j] == dp[i-1][j-1] + match:
                    aln_len[i][j] = aln_len[i-1][j-1] + 1 
            else:
                dp[i][j] = max(dp[i][j], dp[i-1][j-1] + mismatch)
                if dp[i][j] == dp[i-1][j-1] + mismatch:
                    aln_len[i][j] = aln_len[i-1][j-1] + 1
            dp[i][j] = max(dp[i][j], dp[i-1][j] + gap)
            if dp[i][j] == dp[i-1][j] + gap:
                aln_len[i][j] = aln_len[i-1][j] + 1
            dp[i][j] = max(dp[i][j], dp[i][j-1] + gap)
            if dp[i][j] == dp[i][j-1] + gap:
                aln_len[i][j] = aln_len[i][j-1] + 1
            if dp[i][j] == 0:
                aln_len[i][j] = 1
 
    max_i, max_j = 0, 0
    for i in range(1, len(seq2) + 1):
        for j in range(1, len(seq1) + 1):
            if dp[i][j] > dp[max_i][max_j]:
                max_i, max_j = i, j
                # print(i, j, aln_len[i][j], dp[i][j], minlen)
                # alns.append(restore_aln(dp, i, j, match, mismatch, gap, seq1, seq2))
                # print(alns[-1][0])
                # print(alns[-1][1])
                # print("")
    alns.append(restore_aln(dp, max_i, max_j, match, mismatch, gap, seq1, seq2))
    return alns[0], dp[max_i][max_j]


def check_corrupted_as_local(seq, monomers):
    min_ed = len(seq)
    alsat = ""
    max_score = -1
    monomers2 = []
    for m in monomers:
        monomers2.append(m)
        monomers2.append(make_record(m.seq.reverse_complement(), m.id + "'", m.id + "'"))
    for m in monomers2:
        print(m.name)
        aln, score = local_alignments(seq, m.seq)
        print(score)
        print("".join(aln[0]))
        print("".join(aln[1]))
        print("Edist HW,")
        edist_hw([seq, m.seq])
        print("")
        if len("".join(aln[0]).replace("-","")) > max(20, 0.2*len(seq)):
            max_score = max(score, max_score)
    seq = collapse_homo(seq)
    if max_score != -1:
        alsat = "CM"
    else:
        print(len(seq), min_ed)
    return alsat

def restore_aln2(dp, i, j, match, mismatch, gap, s1, s2):
    aln = [[],[]]
    x, y = i, j
    cnt = 0
    while x > 0 or y > 0:
        # print(x, y, dp[x][y], len(s1), len(s2))
        # if x-1 > -1:
        #     print(" ", x-1, y, dp[x-1][y])
        # if y-1 > -1:
        #     print(" ", x, y-1, dp[x][y-1])
        # if x-1 > -1 and y-1 > -1:
        #     print(" ", x-1, y-1)
        #     print(" ", dp[x-1][y-1], s1[y-1], s2[x-1])
        # print("")
        if x-1 > -1 and y -1 > -1 and ((s1[y-1] == s2[x-1] and dp[x][y] == dp[x-1][y-1] + match) or (s1[y-1] != s2[x-1] and dp[x][y] == dp[x-1][y-1] + mismatch)):
            aln[0].append(s1[y-1])
            aln[1].append(s2[x-1])
            x -= 1
            y -= 1
        elif x-1 > -1 and dp[x][y] == dp[x-1][y] + gap:
            aln[0].append("-")
            aln[1].append(s2[x-1])
            x -= 1
        elif y-1 > -1 and dp[x][y] == dp[x][y-1] + gap:
            aln[0].append(s1[y-1])
            aln[1].append("-")
            y -= 1
        cnt += 1
        if cnt > 300:
            print(cnt, len(aln[0]), len(aln[1]))
            exit(-1)
    return [aln[0][::-1], aln[1][::-1]]

def prefix_alignment(seq1, seq2, dp):
    match, mismatch, gap = 0, 1, 1
    alns = []
    dp[0][0] = 0
    for i in range(1, len(seq2) + 1):
        dp[i][0] = dp[i-1][0] + gap
    for i in range(1, len(seq1) + 1):
        dp[0][i] = dp[0][i-1] + gap

    best_scores = []
    w = 5
    for i in range(1, len(seq2) + 1):
        best_score, best_j = 100500, -1
        for j in range(max(1, i - w-1), min(len(seq1) + 1, i + w + 1)):
            dp[i][j] = 100500
        for j in range(max(1, i - w), min(len(seq1) + 1, i + w)):
            if seq1[j - 1] == seq2[i - 1]:
                dp[i][j] = dp[i-1][j-1] + match
            else:
                dp[i][j] = dp[i-1][j-1] + mismatch
            if j < i - 1 + w:
                dp[i][j] = min(dp[i][j], dp[i-1][j] + gap)
            if j - 1 >= i - w:    
                dp[i][j] = min(dp[i][j], dp[i][j-1] + gap)
            if best_score > dp[i][j]:
                best_score, best_j = dp[i][j], j
        best_scores.append([best_score, best_j-1])
    #aln = restore_aln2(dp, max_i, max_j, match, mismatch, gap, seq1, seq2)
    #print("".join(aln[0]))
    #print("".join(aln[1]))    
    return best_scores, dp


def cut_monomers(seqs, monomers, th_mono):
    size = 190
    INF = 100500
    dp = [[INF for _ in range(size)] for _ in range(size)]
    bdp = [[INF for _ in range(size)] for _ in range(size)]
    monomers2 = []
    for m in monomers:
        monomers2.append(m)
        monomers2.append(make_record(m.seq.reverse_complement(), m.id + "'", m.id + "'"))
    for i in range(len(seqs)):
        s = seqs[i]
        if 200 < len(s["seq"]):
            print("i=", i, len(seqs), len(s["seq"]), "read=", s["r"], "start=", s["s"], "end=", s["e"])
            prev_id, prev_seq = s["prev_id"], s["prev_seq"]
            next_id, next_seq = s["next_id"], s["next_seq"]
            cur_seq = s["seq"][:size-1]
            best_ed, best_end1, best_end2, best_end3 = INF, -1, -1, -1
            bcur_seq = s["seq"][-size + 1:]
            best_aln, bbest_aln = None, None
            best_m = None
            for m in monomers2:
                changed = False
                print(len(cur_seq), len(m.seq))
                best_eds, dp = prefix_alignment(cur_seq, m.seq, dp)
                bbest_eds, bdp = prefix_alignment(bcur_seq[::-1], m.seq[::-1], bdp)
                for k in range(len(best_eds)):
                    if best_eds[k][0] + bbest_eds[len(m.seq) - k - 1][0] < best_ed:
                        changed = True
                        best_ed = best_eds[k][0] + bbest_eds[len(m.seq) - k - 1][0]
                        best_end1, best_end2, best_end3 = k, best_eds[k][1], len(s["seq"]) - bbest_eds[len(m.seq) - k - 1][1]

                if changed:
                    best_m = m.name
                    best_aln = restore_aln2(dp, best_end1 + 1, best_end2 + 1, 0, 1, 1, cur_seq, m.seq)
                    bbest_aln = restore_aln2(bdp, len(m.seq) - best_end1, len(s["seq"]) - best_end3 + 1, 0, 1, 1, bcur_seq[::-1], m.seq[::-1])
            print(best_m)
            print(best_ed, (171 - best_ed)/171*100, best_end1, best_end2, best_end3, th_mono)
            print("".join(best_aln[0]))
            print("".join(best_aln[1]))
            print("")
            print("".join(bbest_aln[0]))
            print("".join(bbest_aln[1]))    
            print("")
            print(len(s["seq"]), len(seqs[i]["seq"]))
            if (171 - best_ed)/171*100 > th_mono:
                seqs[i]["seq"] = s["seq"][best_end2: best_end3]
            print(len(seqs[i]["seq"]))
            print()
        # else:
        #     seqs[i]["seq"] = s["seq"]

    #exit(-1)
    return seqs

def find_set(ind, parent):
    if ind == parent[ind]:
        return ind
    parent[ind] = find_set(parent[ind], parent)
    return parent[ind]

def union_sets(a, b, parent, rank):
    a = find_set(a, parent)
    b = find_set(b, parent)
    if a != b:
        if rank[a] < rank[b]:
            a, b = b, a
        parent[b] = a
        if rank[a] == rank[b]:
            rank[a] += 1

def collapse_homo(sequence):
    new_sequence = sequence[0]
    for i in range(1, len(sequence)):
        if sequence[i] != sequence[i-1]:
            new_sequence += sequence[i]
    return Seq(new_sequence)

def cluster_by_outer_ed(sequences, th=20):
    clusters = []
    parent = [i for i in range(len(sequences))]
    rank = [0 for _ in range(len(sequences))]
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            #h_sequence_i, h_sequence_j = collapse_homo(sequences[i]), collapse_homo(sequences[j])
            h_sequence_i, h_sequence_j = sequences[i]["seq"], sequences[j]["seq"]
            ed = edist([h_sequence_i, h_sequence_j])
            min_ed = min(ed, edist([h_sequence_i.reverse_complement(), h_sequence_j]))
            if min_ed*100/min(len(h_sequence_j), len(h_sequence_i)) <= th:
                union_sets(i, j, parent, rank)

    clusters_id = {}
    for i in range(len(sequences)):
        ind = find_set(i, parent)
        ed = edist([sequences[i]["seq"], sequences[ind]["seq"]])
        min_ed = min(ed, edist([sequences[i]["seq"].reverse_complement(), sequences[ind]["seq"]]))
        if ind not in clusters_id:
            clusters_id[ind] = []
        if min_ed != ed:
            sequences[i]["rev"] = True
            clusters_id[ind].append(sequences[i])
        else:
            clusters_id[ind].append(sequences[i])

    for cl in clusters_id:
        clusters.append(clusters_id[cl])
        if len(clusters[-1]) > 1:
            print(len(clusters), len(clusters[-1]), len(clusters[-1][0]["seq"]), clusters[-1][0]["r"])
    return clusters

def construct_representative(clusters, reads, monomers, minlen):
    res = []
    new_clusters = []
    for i in range(len(clusters)):
        cl = clusters[i]
        best_ind, best_score = 0, len(cl[0]["seq"])
        for j in range(len(cl)):
            cur_score = -1
            for k in range(len(cl)):
                if j != k:
                    a, b = cl[k]["seq"], cl[j]["seq"]
                    if cl[k]["rev"] != cl[j]["rev"]:
                        b = b.reverse_complement()
                    ed = edist([a, b])
                    if cur_score == -1 or cur_score < ed:
                        cur_score = ed
            if cur_score < best_score:
                best_score, best_ind = cur_score, j
        if len(cl) > 1:
            print(i + 1, best_ind, best_score, cl[best_ind]["r"], cl[best_ind]["s"], cl[best_ind]["e"])

        rev = cl[best_ind]["rev"]
        alsat = ""
        if len(cl[best_ind]["seq"]) <= 200:
            print("Checking... ", "NM_" + str(i + 1) + "_" + str(len(cl)) + "_"  + str(len(cl[best_ind]["seq"])))
            alsat = check_corrupted_as(cl[best_ind]["seq"], monomers)
        if len(alsat) == 0:
            name = "NM_" + str(i + 1) + "_" + str(len(cl)) + "_"  + str(len(cl[best_ind]["seq"]))
        else:
            name = "CM_" + str(i + 1) + "_" + str(len(cl)) + "_"  + str(len(cl[best_ind]["seq"]))
        res.append(make_record(cl[best_ind]["seq"], name, name))
        if rev:
            for j in range(len(cl)):
                clusters[i][j]["rev"] = not clusters[i][j]["rev"]
    return res, clusters


def contruct_nm_regions(decomposition, reads, monomers, params):
    reads_regions = load_regions(decomposition)
    nonmono_seq, clusters = identify_nm(reads_regions, params, reads, monomers)
    return nonmono_seq, clusters


def identify_nm(reads_regions, params, reads, monomers, datatype):
    regions = extract_nonmono_regions(reads_regions, reads, datatype, params["th_border"], params["th_mono"])
    cnt = 0
    seqs = []
    reads_with_regions = set()
    for r in regions:
        for region in regions[r]:
            if 20000 > region["e"] -  region["s"] + 1 > params["min_len"]:
                #print(r, len(regions[r]), region["s"],  region["e"], region["e"] -  region["s"])
                #print(region["seq"])
                seqs.append({"seq": Seq(region["seq"]), "r": r, "s": region["s"], "e": region["e"], "rev": False\
                            ,"prev_id": region["prev_id"], "prev_seq": region["prev_seq"], "next_id": region["next_id"], "next_seq": region["next_seq"]})
                cnt += 1
                reads_with_regions.add(r)

    print("Found", cnt, "potential non-monomeric regions in", len(reads_with_regions), "reads")
    print("Removing monomer ends..")
    seqs = cut_monomers(seqs, monomers, params["th_mono"])
    print("Clustering regions..")
    clusters = cluster_by_outer_ed(seqs, params["ed"])
    print("Identify representatives..")
    consensus, clusters = construct_representative(clusters, reads, monomers, params["min_len"])
    total = len(consensus)
    consensus = sorted(consensus, key = lambda x: -int(x.id.split("_")[3]))
    if datatype != "asm":
        consensus = [x for x in consensus if int(x.id.split("_")[2]) > 1]
        filtered = len(consensus)
        print("Number of non-monomeric elements ", total, ". With more than one occurence ", filtered)
    else:
        print("Number of non-monomeric elements ", total)
    # for c in consensus:
    #     print(">" + c.id)
    #     print(c.seq)
    # for i in range(len(consensus)):
    #     for j in range(i + 1, len(consensus)):
    #         ed = edist_hw([consensus[i].seq, consensus[j].seq])
    #         min_ed = min(ed, edist_hw([consensus[i].seq.reverse_complement(), consensus[j].seq]))
    #         min_ed = min(min_ed, edist_hw([consensus[j].seq.reverse_complement(), consensus[i].seq]))
    #         min_ed = min(min_ed, edist_hw([consensus[j], consensus[i].seq]))
    #         print(consensus[i].id, consensus[j].id, min_ed)
    return consensus, clusters


def form_nm_decomposition(non_mono, clusters, reads, decomposition_file, outfile):
    non_mono_byreads = {}
    for m in non_mono:
        cl_id = int(m.id.split("_")[1])
        for region in clusters[cl_id - 1]:
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
    with open(outfile, "w") as fout:
        with open(decomposition_file, "r") as fin:
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Searches for non-monomeric regions')
    parser.add_argument('sequences', help='fasta-file with long reads or genomic sequences')
    parser.add_argument('monomers', help='fasta-file with monomers')
    parser.add_argument('decomposition', help='tsv-file with sequences decomposition')
    parser.add_argument('output', help='tsv-file to save new decomposition')
    parser.add_argument('-d', '--data-type',  help='type of reads (hifi or ont)', choices=["hifi", "ont", "asm"], required = True)
    parser.add_argument('--min-monomer',  help='minimum identity of monomer (70 for ONT and 85 for Hifi, by default)', type=int, default=-1, required = False)
    parser.add_argument('--min-reliable',  help='minimum identity of reliable monomer (95, by default)', type=int, default=95, required = False)
    parser.add_argument('--min-length',  help='minimum length of non-monomeric region (200 bp, by default)', type=int, default=200, required = False)
    parser.add_argument('--max-diff',  help='threshold on identiy for clustering (20 for ONT and 5 for Hifi, by default)', type=int, default=-1, required = False)
    #parser.add_argument('--use-clustal',  help='use clustal to construct representatives', action="store_true")
    args = parser.parse_args()

    params = {"th_border": args.min_reliable, "th_mono": args.min_monomer, "min_len": args.min_length, "ed": args.max_diff}
    if params["ed"] == -1:
        if args.data_type == "hifi" or args.data_type == "asm":
            params["ed"] = 5
        else:
            params["ed"] = 20

    if params["th_mono"] == -1:
        if args.data_type == "hifi" or args.data_type == "asm":
            params["th_mono"] = 75
        else:
            params["th_mono"] = 70

    reads = load_fasta(args.sequences, "map")
    monomers = load_fasta(args.monomers)
    all_regions = load_regions(args.decomposition)
    non_mono, clusters = identify_nm(all_regions, params, reads, monomers, args.data_type)
    new_monomer_file = args.output[:-len(".tsv")] + "_monomers_with_nm_regions.fasta"
    new_dec_elements = monomers + non_mono
    print("Saving new set of elements to decompose to ", new_monomer_file)
    save_fasta(new_monomer_file, new_dec_elements)
    print("Saving new decomposition to ", args.output)
    form_nm_decomposition(non_mono, clusters, reads, args.decomposition, args.output)
