import sys
import edlib

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline

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
    result = edlib.align(str(lst[0]), str(lst[1]), mode="HW")
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
            sseqid = sseqid.split()[0]
            if sseqid not in reads_regions:
                reads_regions[sseqid] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            reads_regions[sseqid].append({"idnt": idnt, "s": s, "e": e})
    for r in reads_regions:
        reads_regions[r] = sorted(reads_regions[r], key=lambda x: x["s"])
    return reads_regions

def extract_nonmono_regions(reads_regions, reads, th_border = 95, th_mono = 85):
    window = 2
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
            for i in range(len(reads_regions[r])):
                if reads_regions[r][i]["idnt"] < th_mono:
                    if in_region:
                        collapsed_nonmono_regions[r][-1].append(reads_regions[r][i])
                    else:
                        in_region = True
                        #print("new region")
                        collapsed_nonmono_regions[r].append([reads_regions[r][i]])
                    #print(reads_regions[r][i]["s"], reads_regions[r][i]["e"],reads_regions[r][i]["idnt"], th_mono)
                else:
                    in_region = False
    new_collapsed_nonmono_regions = {}
    for r in collapsed_nonmono_regions:
        if len(collapsed_nonmono_regions[r]) > 0:
            new_collapsed_nonmono_regions[r] = []
            for region in collapsed_nonmono_regions[r]:
                new_collapsed_nonmono_regions[r].append({ "seq": "".join(reads[r].seq[region[0]["s"] : region[-1]["e"] + 1]), "s": region[0]["s"], "e": region[-1]["e"] })
    return new_collapsed_nonmono_regions

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

def construct_representative(clusters, reads):
    res = []
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
        name = "NM_" + str(i + 1) + "_" + str(len(cl)) + "_"  + str(len(cl[best_ind]["seq"]))
        res.append(make_record(cl[best_ind]["seq"], name, name))
        if rev:
            for j in range(len(cl)):
                clusters[i][j]["rev"] = not clusters[i][j]["rev"]
    return res, clusters

def save_clusters(clusters, prefix):
    names = []
    for i in range(len(clusters)):
        sequences = []
        for j in range(len(clusters[i])):
            sequences.append(make_record(clusters[i][j], "cl" + str(i) + "it" + str(j), "cl" + str(i) + "it" + str(j)))
        num = len(sequences)
        if len(sequences) == 1:
            sequences.append(sequences[0])
        save_fasta(prefix + "/cl" +  str(i) + "_" + str(num) + ".fasta", sequences)
        names.append(prefix + "/cl" +  str(i) + "_" + str(num) + ".fasta")
    return names

def construct_consensus(clusters):
    clustalw_exe = r"/home/tdvorkina/soft/clustalo-1.2.4-Ubuntu-x86_64"
    for cl in clusters:
        cline = ClustalwCommandline(clustalw_exe, infile=cl, outfile=cl[:-len("fasta")] + "clu")
        stdout, stderr = cline()

    ind = 0
    seq_consensus = []
    for cl in clusters:
        filename = cl[:-len("fasta")] + "clu"
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
        s_consensus = ""
        for i in range(len(total_alns[0])):
            score = {"A": 0, "C": 0, "G": 0, "T": 0, "-": 0}
            for j in range(len(total_alns)):
                score[total_alns[j][i]] += 1
            max_score = "A"
            for s in score:
                if score[s] > score[max_score]:
                    max_score = s
            if max_score != "-":
                s_consensus += max_score
        #print(cl, s_consensus)
        name = cl.split("/")[-1].split(".")[0]
        seq_consensus.append(make_record(Seq(s_consensus), name + "_" + str(ind)  + "_" + str(len(cl))  + "_" + str(len(s_consensus)), name + "_" + str(ind)  + "_" + str(len(cl))  + "_" + str(len(s_consensus))))
        i += 1
    return seq_consensus


def contruct_nm_regions(decomposition, reads, params):
    reads_regions = load_regions(decomposition)
    nonmono_seq, clusters = identify_nm(reads_regions, params, reads)
    return nonmono_seq, clusters


def identify_nm(reads_regions, params, reads):
    regions = extract_nonmono_regions(reads_regions, reads, params["th_border"], params["th_mono"])
    cnt = 0
    seqs = []
    reads_with_regions = set()
    for r in regions:
        for region in regions[r]:
            if 20000 > region["e"] -  region["s"] > params["min_len"]:
                #print(r, len(regions[r]), region["s"],  region["e"], region["e"] -  region["s"])
                #print(region["seq"])
                seqs.append({"seq": Seq(region["seq"]), "r": r, "s": region["s"], "e": region["e"], "rev": False})
                cnt += 1
                reads_with_regions.add(r)

    print("Found", cnt, "potential non-monomeric regions in", len(reads_with_regions), "reads")
    print("Clustering regions..")
    clusters = cluster_by_outer_ed(seqs, params["ed"])
    print("Identify representatives..")
    consensus, clusters = construct_representative(clusters, reads)
    total = len(consensus)
    consensus = sorted(consensus, key = lambda x: -int(x.id.split("_")[3]))
    consensus = [x for x in consensus if int(x.id.split("_")[2]) > 1]
    filtered = len(consensus)
    print("Number of non-monomeric elements ", total, ". With more than one occurence ", filtered)
    for c in consensus:
        print(">" + c.id)
        print(c.seq)
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
    parser.add_argument('-d', '--data-type',  help='type of reads (hifi or ont)', choices=["hifi", "ont"], required = True)
    parser.add_argument('--min-monomer',  help='minimum identity of monomer (70 for ONT and 85 for Hifi, by default)', type=int, default=-1, required = False)
    parser.add_argument('--min-reliable',  help='minimum identity of reliable monomer (95, by default)', type=int, default=95, required = False)
    parser.add_argument('--min-length',  help='minimum length of non-monomeric region (200 bp, by default)', type=int, default=200, required = False)
    parser.add_argument('--max-diff',  help='threshold on identiy for clustering (20 for ONT and 5 for Hifi, by default)', type=int, default=-1, required = False)
    #parser.add_argument('--use-clustal',  help='use clustal to construct representatives', action="store_true")
    args = parser.parse_args()

    params = {"th_border": args.min_reliable, "th_mono": args.min_monomer, "min_len": args.min_length, "ed": args.max_diff}
    if params["ed"] == -1:
        if args.data_type == "hifi":
            params["ed"] = 5
        else:
            params["ed"] = 20

    if params["th_mono"] == -1:
        if args.data_type == "hifi":
            params["th_mono"] = 85
        else:
            params["th_mono"] = 70

    reads = load_fasta(args.sequences, "map")
    all_regions = load_regions(args.decomposition)
    non_mono, clusters = identify_nm(all_regions, params, reads)
    # if args.use_clustal:
    #     print("Saving clusters to", prefix)
    #     filenames = save_clusters(clusters, prefix)
    #     print("Constructing consensus using clustal..")
    #     consensus = construct_consensus(filenames)
    # else:
    monomers = load_fasta(args.monomers)
    new_monomer_file = args.output[:-len(".tsv")] + "_monomers_with_nm_regions.fasta"
    new_dec_elements = monomers + non_mono
    print("Saving new set of elements to decompose to ", new_monomer_file)
    save_fasta(new_monomer_file, new_dec_elements)
    print("Saving new decomposition to ", args.output)
    form_nm_decomposition(non_mono, clusters, reads, args.decomposition, args.output)
