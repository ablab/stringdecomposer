#!/usr/bin/env python
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import pysam
import edlib

import sys, os

def cnt_identity(lst, cur_mode = "NW"):
    if len(str(lst[0])) == 0:
        return 0
    if len(str(lst[1])) == 0:
        return 0
    result = edlib.align(str(lst[0]), str(lst[1]), mode=cur_mode, task="distance")
    return 100 - result["editDistance"]*100//max(len(str(lst[0])), len(str(lst[1])))

def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    else:
        records = list(SeqIO.parse(filename, "fasta"))
    return records

def load_blast_tsv(filename, identity):
    reads_mapping = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            qseqid, sseqid, pident, qlen, slen, length, qstart, qend, sstart, send = ln.strip().split("\t")
            if int(qlen)*0.8 < int(qend) - int(qstart) + 1 and float(pident) > identity:
                if sseqid not in reads_mapping:
                    reads_mapping[sseqid] = []
                s, e = int(sstart), int(send)
                rev = False
                if s > e:
                    s, e = e, s
                    rev = True
                reads_mapping[sseqid].append({"qid": qseqid, "s": s, "e": e, "rev":rev})
    for r in reads_mapping:
        reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: x["s"])
    return reads_mapping

def load_minimap2_tsv(filename, reads, monomers, identity):
    reads_mapping = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            qseqid, _, qstart, qend, st, sseqid, _, sstart, send  = ln.strip().split("\t")[:9]
            if int(send) - int(sstart) > 0.8*len(monomers[sseqid].seq):
                qseq = reads[qseqid].seq[int(qstart):int(qend)]
                if st == '-':
                    qseq = qseq.reverse_complement()
                if cnt_identity([qseq, monomers[sseqid].seq]) > identity:
                    # print str(qseq)
                    # print str(monomers[sseqid].seq)
                    # print ""
                    if qseqid not in reads_mapping:
                        reads_mapping[qseqid] = []
                    s, e = int(qstart), int(qend)
                    rev = False
                    if s > e:
                        s, e = e, s
                        rev = True
                    reads_mapping[qseqid].append({"qid": sseqid, "s": s, "e": e, "rev": rev})
    for r in reads_mapping:
        reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: x["s"])
    return reads_mapping


def load_decomposition_tsv(filename, reads, monomers):
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
            if q == "+":
                if rev:
                    idnt = cnt_identity([reads[sseqid].seq[s: e + 1], monomers[qseqid].seq.reverse_complement()])
                else:
                    idnt = cnt_identity([reads[sseqid].seq[s: e + 1], monomers[qseqid].seq])
                reads_mapping[sseqid].append({"qid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt})
    for r in reads_mapping:
        if len(reads_mapping[r]) > 0 and reads_mapping[r][0]["rev"]:
            reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (-x["e"], x["s"]))
        else:
            reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (x["s"], -x["e"]))
    new_mapping = {}
    for r in reads_mapping:
        new_mapping[r] = []
        if len(reads_mapping[r]) > 0:
            if reads_mapping[r][0]["rev"]:
                new_mapping[r].append(reads_mapping[r][0])
                for i in range(len(reads_mapping[r]) - 1):
                    if not (reads_mapping[r][i]["e"] >= reads_mapping[r][i + 1]["e"] - 5 and reads_mapping[r][i]["s"] - 5 <= reads_mapping[r][i + 1]["s"]):
                        new_mapping[r].append(reads_mapping[r][i + 1])
                        # if reads_mapping[r][i]["qid"] == reads_mapping[r][i + 1]["qid"]:
                        #     print([reads_mapping[r][i], reads_mapping[r][i+1]]) 
            else:
                new_mapping[r].append(reads_mapping[r][0])
                for i in range(len(reads_mapping[r]) - 1):
                    if not (reads_mapping[r][i]["s"] -5 <= reads_mapping[r][i + 1]["s"] and reads_mapping[r][i]["e"] >= reads_mapping[r][i + 1]["e"] - 5):
                        new_mapping[r].append(reads_mapping[r][i + 1])
                        # if reads_mapping[r][i]["qid"] == reads_mapping[r][i + 1]["qid"]:
                        #     print([reads_mapping[r][i], reads_mapping[r][i+1]])

    return new_mapping

def load_assembly_reads():
    filename = "assembly_reads.txt"
    names = set()
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            names.add(ln.strip())
    return names
    

def analyze_order(read_hits, monomers, reads=[]):
    order_lst = {}
    for m in monomers:
        order_lst[m] = {} 
        for m2 in monomers:
            order_lst[m][m2] = 0
    bad_reads = 0
    for r in sorted(read_hits.keys()):
        if len(reads) == 0 or r in reads:
            good_read = True
            rev = False
            for i in range(len(read_hits[r]) - 1):
                h1 = read_hits[r][i]
                h2 = read_hits[r][i + 1]
                rev = h1["rev"]
                if h1["rev"] != h2["rev"]:
                    good_read = False
                    break
            if good_read:  
                for i in range(len(read_hits[r]) - 1):
                    h1 = read_hits[r][i]
                    h2 = read_hits[r][i + 1]
                    if  (not rev and abs(h2["s"] - h1["e"]) < 10) or (rev and abs(h1["s"] - h2["e"]) < 10):
                        order_lst[h1["qid"]][h2["qid"]] += 1
                        #s1, s2 = int(h1["qid"].split("-")[1]), int(h2["qid"].split(":")[1].split("-")[0])
                        s1, s2 = int(h1["qid"].split("_")[1]), int(h2["qid"].split("_")[1]) 
                        if s1 + 1 != s2 and not (s1 == 12 and s2 == 1):
                            print(r + " " + h1["qid"] + " " + str(h1["s"]) + " " + str(h1["e"]) + " " + h2["qid"] + " " + str(h2["s"]) + " " + str(h2["e"]))
                    else:
                        if rev:
                            print(r + " " + h1["qid"] + " " + str(h1["s"]) + " " + str(h1["e"]) + " " + h2["qid"] + " " + str(h2["s"]) + " " + str(h2["e"]) + " " + str(h1["s"] - h2["e"]))
                        else:
                            print(r + " " + h1["qid"] + " " + str(h1["s"]) + " " + str(h1["e"]) + " " + h2["qid"] + " " + str(h2["s"]) + " " + str(h2["e"]) + " " + str(h2["s"] - h1["e"]))
            else:
                bad_reads += 1

    #print("Bad reads " + str(bad_reads))

    symbol = {}
    for m in monomers:
        #symbol[m] = int(m.split(":")[1].split("-")[0])
        symbol[m] = int(m.split("_")[1])

    s_r = {}
    for s in symbol:
        s_r[symbol[s]] = s

    for m in sorted([symbol[k] for k in order_lst.keys()]):
        ans = []
        for m2 in sorted([symbol[kk] for kk in order_lst[s_r[m]].keys()]):
            ans.append(str(order_lst[s_r[m]][s_r[m2]]))
        print("\t".join(ans))

def identity(read_hits, reads, monomers):
    sum_identity = 0.0
    cnt = 0
    for r in read_hits:
        for h in read_hits[r]:
            cnt += 1
            identity = max(cnt_identity([monomers[h["qid"]].seq, reads[r].seq[h["s"]: h["e"]] ]), \
                        cnt_identity([monomers[h["qid"]].seq, reads[r].seq[h["s"]: h["e"]].reverse_complement() ]))
            sum_identity += identity
    print("Avg identity =" + str(sum_identity//cnt))

def hor_coverage(read_hits, reads):
    num = 0
    for r in read_hits:
        # good_read = True
        # for i in range(len(read_hits[r]) - 1):
        #     h1 = read_hits[r][i]
        #     h2 = read_hits[r][i + 1]
        #     if h1["rev"] != h2["rev"]:
        #         good_read = False
        #         break
        # if good_read: 
        if len(read_hits[r]) >= 36:
            num += 1
    print(num)

def analyze_ncrf(filename):
    total_sz, aln_sz, sum_identity, cnt = 0, 0, 0, 0
    n = 0
    reads = {}
    with open(filename, "r") as fin:
        read_name = None
        for ln in fin.readlines():
            lst = ln.strip().split()
            if len(lst) > 4:
                #print(lst[:4])
                read_name, rlen, alen, c1, c2, aln_seq1 = lst[0], int(lst[1]), int(lst[2][:-2]), int(lst[3].split("-")[0]), int(lst[3].split("-")[1]), lst[4]
                aln_seq1 = aln_seq1.upper().replace("-", "")
                reads[read_name] = rlen
                #print(read_name)
            elif len(lst) > 3:
                #print(lst[:3])
                name, rplen, aln_seq2 = lst[0], int(lst[1][:-2]), lst[3]
                aln_seq2 = aln_seq2.upper().replace("-", "")
                
                if name.endswith("-"):
                    aln_seq2 = str(Seq(aln_seq2).reverse_complement())
                aln_sz += alen
                n += 1
                #print([len(aln_seq1[c1:c2]), len(aln_seq2)])
                idd = cnt_identity([aln_seq1, aln_seq2])
                #if idd < 50:
                    # print(len(aln_seq1))
                    # print(aln_seq1[c1:c2])
                    # print([read_name, c1, c2, idd])
                    #break
                sum_identity += idd
                cnt += 1
    for r in reads:
        total_sz += reads[r]
    print("Num reads " + str(len(reads)) + " len: " + str(total_sz))   
    print("Percent of covered: " + str(aln_sz*100//total_sz) + "%" + " Total: " + str(total_sz))
    print("Average identity: " + str(sum_identity//cnt))


def coverage(read_hits, reads):
    props = []
    sum_len = 0
    sum_filled = 0
    reads_num = 0
    bad_reads = 0
    bad_len = 0
    avg_idnt, mapping_num = 0, 0
    for r in read_hits:
        good_read = True
        for i in range(len(read_hits[r]) - 1):
            h1 = read_hits[r][i]
            h2 = read_hits[r][i + 1]
            if h1["rev"] != h2["rev"]:
                good_read = False
                break
        if good_read:  
            read_len = len(reads[r].seq)
            mappings = read_hits[r]
            ar = [0 for _ in range(read_len)]
            for m in mappings:
                s, e = m["s"], m["e"]
                for i in range(s-1, e):
                    ar[i] = 1
                avg_idnt += m["idnt"]
                mapping_num += 1
            cnt = sum(ar)
            props.append(cnt*100//read_len)
            sum_len += read_len
            sum_filled += cnt
            reads_num += 1
        else:
            bad_reads += 1
            sum_len += len(reads[r].seq)
            bad_len += len(reads[r].seq)
    print("Bad reads " + str(bad_reads))
    for r in reads:
        if r not in read_hits:
            print(r)
    print("Bad reads all " + str(bad_reads) + " " + str(bad_len) + " " + str(bad_len*100//sum_len))
    sprops = sorted(props)
    # print(sprops[:10])
    # print(sprops[-10:])
    # print(sprops[len(sprops)/2])
    print("Percent of covered " + str(sum_filled*100//sum_len) + "%") #34%
    print("Avg identity " + str(avg_idnt//mapping_num) + "%")
    print("Filled " + str(sum_filled) + " Total length " + str(sum_len))
    print("Covered reads " + str(reads_num) + " Total read num " + str(len(reads)))


#analyze_ncrf("/home/tdvorkina/projects/centroFlye/cenX_v0_8_2/report_reads.ncrf")
#analyze_ncrf("/home/tdvorkina/projects/centroFlye/cenX_v0_8_2/report_assembly.ncrf")
#analyze_ncrf("/home/tdvorkina/projects/centroFlye/chr6/report.ncrf")

#reads = load_fasta("./cenX_centroFlye_assembly_v0_8_1.fasta", "map")
#monomers = load_fasta("../chr6/D6Z1_hg38_monomers.fasta", "map")

# minimap2_hits = load_minimap2_tsv("./minimap/aln_ref.paf", reads, monomers, 80)
# blast_hits = load_blast_tsv("blast_assembly/monomers_mapping.tsv", 80)
#dec_hits = load_decomposition_tsv("decomposition_ref_test2.tsv")
# print("\nref: Minimap2 hits")
# identity(minimap2_hits, reads, monomers)
# coverage(minimap2_hits, reads)
# print("\nref: BLAST hits")
# identity(dec_hits, reads, monomers)
# coverage(blast_hits, reads)
# print("\nref: Decomposition hits")
# identity(dec_hits, reads, monomers)
# coverage(dec_hits, reads)


#reads = load_fasta("../chr6/chm13_d6z1.rds_ge100kb.NHGRI.guppy_3.1.5_fixed.fa", "map")
#blast_hits = load_blast_tsv("blast/monomers_mapping_50.tsv", 80)
#minimap2_hits = load_minimap2_tsv("./minimap/aln_reads.paf", reads, monomers, 80)
#dec_hits = load_decomposition_tsv("../chr6/decomposition_reads_new.tsv")
#coverage(dec_hits, reads)
#reads = load_assembly_reads()
#analyze_order(dec_hits, monomers)

# print("\nreads: Minimap2 hits")
# #identity(minimap2_hits, reads, monomers)
# hor_coverage(minimap2_hits, reads)
# print("\nreads: BLAST hits")
# #identity(minimap2_hits, reads, monomers)
# hor_coverage(blast_hits, reads)
# print("\nreads: Decomposition hits")
# #identity(minimap2_hits, reads, monomers)
# hor_coverage(dec_hits, reads)



reads = load_fasta("../chrX/results/centromeric_reads/centromeric_reads.fasta", "map")
monomers = load_fasta("../chrX/DXZ1_rc_star_monomers.fasta", "map")
#dec_hits = load_decomposition_tsv("../chrX/decomposition_reads_star_new.tsv", reads, monomers)
dec_hits = load_decomposition_tsv("../chrX/decomposition_reads_star_new.tsv", reads, monomers)
coverage(dec_hits, reads)
analyze_order(dec_hits, monomers)