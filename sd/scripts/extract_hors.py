from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import sys
import os
from os import listdir
from os.path import isfile, isdir, join
import argparse

import logging

import re

sd_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.append(sd_path)
from sd.utils.bio import read_bio_seqs

logger = logging.getLogger("SD.scripts.extract_hors")

def load_dec(filename, min_idnt, min_reliable):
    reads_mapping = {}
    ss = ""
    monomers = set()
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt  = ln.strip().split("\t")[:5]
            sseqid = sseqid.split()[0]
            if sseqid not in reads_mapping:
                    reads_mapping[sseqid] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            rev = False
            if qseqid.endswith("'"):
                rev = True
                qseqid = qseqid[:-1]
            qseqid = qseqid.split()[0]
            qseqid = qseqid.split("(")[0]
            if idnt >= min_idnt:
                monomers.add(qseqid)
                reads_mapping[sseqid].append({"qid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt})
            else:
                reads_mapping[sseqid].append({"qid": "NM","s": s, "e": e, "rev": False, "idnt": idnt})

    new_reads_mapping = {}
    cnt = 0
    for r in reads_mapping:
        cur_mapping = []
        inside_good = False
        left, right = -1, -1
        for i in range(len(reads_mapping[r])):
            if left == -1 and reads_mapping[r][i]["idnt"] > min_reliable:
                left = i
            if reads_mapping[r][i]["idnt"] > min_reliable:
                right = i
        if len(reads_mapping[r][left: right + 1]) > 36:
            new_reads_mapping[r] = sorted(reads_mapping[r][left: right + 1], key = lambda x: x["s"])
            cnt += right -left+ 1
    return new_reads_mapping


def convert_dec_to_internal_monomers(reads_mapping, monomers):
    monomers_lst = [x for x in monomers]
    monomers_mp = {}
    monomers_mp_r = {}
    print("Monomer mapping: ")
    for i in range(len(monomers_lst)):
        qid = monomers_lst[i]
        if qid != "NM":
            new_qid = "m" + str(i + 1)
        else:
            new_qid = "f" + str(i + 1)
        print(" ", "\t".join([qid, new_qid]))
        monomers_mp[qid] = new_qid
        monomers_mp_r[new_qid] = qid
        monomers_mp_r[new_qid + "'"] = qid + "'"
    monomers_mp["NM"] = "NM"
    monomers_mp_r["NM"] = "NM"
    for c in reads_mapping:
        for j in range(len(reads_mapping[c])):
            if reads_mapping[c][j]["rev"]:
                reads_mapping[c][j]["qid"] = monomers_mp[reads_mapping[c][j]["qid"]] + "'"
            else:
                reads_mapping[c][j]["qid"] = monomers_mp[reads_mapping[c][j]["qid"]]
    
    return reads_mapping, monomers_mp, monomers_mp_r

def build_known_hors(filename, monomers_mp):
    known_hors = []
    known_hors_initial = set()
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            known_hors_initial.add(ln.strip())
            lst = ln.strip().split(",")
            kh = []
            for c in lst:
                kh.append(monomers_mp[c] + "[1]")
            known_hors.append("_".join(kh))
    return known_hors, known_hors_initial

def build_full_hor(new_hor, hors, name):
    new_hor_represent = []
    for c in new_hor.split("_"):
        cc, cnt = c.split("[")[0], int(c.split("[")[1][:-1])
        for j in range(cnt):
            if cc in hors:
                new_hor_represent.append(hors[cc])
            else:
                new_hor_represent.append(cc + "[1]")
    hors[name] = "_".join(new_hor_represent)
    return hors

def known_hors_annotation(reads_annotation, known_hors, hors, hors_lst, h_cnt, min_cnt):
    hors_log = []
    rev_comp_hors = []
    for kh in known_hors:
        rev_comp_hors.append("_".join([x[:-len("[1]")] + "'[1]" for x in kh.split("_")][::-1]))
    for kh in known_hors + rev_comp_hors:
        cnt = 0
        set_size, new_set_size = 0, 0
        for r in reads_annotation:
            annotation = reads_annotation[r]
            annotation_seq = []
            for a in annotation:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            annotation_str = "_".join(annotation_seq)
            set_size += len(annotation_seq)
            annotation_new_lst = annotation_str.split(kh)
            annotation_new_str = annotation_str.replace(kh, "X")
            annotation_new_str = annotation_new_str.replace("X_X", "X")
            new_set_size += len(annotation_new_str.split("_"))
            cnt += len(annotation_new_lst) - 1
        if cnt == 0:
            continue
        h_cnt += 1
        name = "h" + str(h_cnt)
        hors= build_full_hor(kh, hors, name)
        hors_lst.append([name, kh])
        hors_log.append([name, kh.replace("[1]",""), hors[name].replace("[1]",""), str(len(hors[name].split("_"))), str(cnt), \
                                                     str(set_size - new_set_size), str(new_set_size) ])
        print("\t".join([name, kh.replace("[1]",""), hors[name].replace("[1]",""), str(len(hors[name].split("_"))), str(cnt), \
                                                     str(set_size - new_set_size), str(new_set_size) ]), flush=True)
        for r in reads_annotation:
            annotation = reads_annotation[r]
            annotation_seq = []
            for a in annotation:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            annotation = update_annotation(annotation, annotation_seq, kh, name)
            reads_annotation[r] = collapse_annotation(annotation)

    return reads_annotation, hors, hors_lst, h_cnt, hors_log

def update_annotation(annotation, annotation_seq, new_hor, name):
    new_annotation = []
    k = len(new_hor.split("_"))
    i = 0
    while i < len(annotation):
        if i + k < len(annotation):
            cur_seq = "_".join(annotation_seq[i:i+k])
        else:
            cur_seq = ""
        if cur_seq == new_hor:
            new_annotation.append([name, 1, {"s": [annotation[i][2]["s"][0]], "e": [annotation[i + k - 1][2]["e"][-1]]} ])
            i += k
        else:
            new_annotation.append(annotation[i])
            i += 1
    return new_annotation

def collapse_annotation(annotation):
    new_annotation = []
    prev, cnt, start, end = "", 0, [], []
    for i in range(len(annotation)):
        c, c_cnt = annotation[i][:2]
        if c == prev:
            cnt += c_cnt
            start.extend(annotation[i][2]["s"])
            end.extend(annotation[i][2]["e"])
        else:
            if prev != "":
                new_annotation.append([prev, cnt, {"s": start, "e": end}])
            prev = c 
            cnt = c_cnt
            start = annotation[i][2]["s"]
            end = annotation[i][2]["e"]
    i = len(annotation) - 1
    new_annotation.append([prev, cnt, {"s": start, "e": end}])
    return new_annotation

def find_potential_hors(annotation, annotation_seq, min_hor_len, max_hor_len, hors, potential_hors_all, potential_hors_names_all, set_size):
    potential_hors = {}
    potential_hors_names = []
    annotation_str = "_".join(annotation_seq)
    for i in range(len(annotation)):
        end_ind = i
        len_subseq, subseq = 0, []
        len_monomer_subseq = 0
        while end_ind < len(annotation) and len_subseq < max_hor_len \
              and annotation[end_ind][0] != "NM" and not annotation[end_ind][0].startswith("f") :
            len_subseq += annotation[end_ind][1]
            if annotation[end_ind][0] in hors:
                len_monomer_subseq += len(hors[annotation[end_ind][0]].split("_"))*annotation[end_ind][1]
            else:
                len_monomer_subseq += annotation[end_ind][1]
            subseq.append(annotation_seq[end_ind])
            end_ind += 1
            if min_hor_len < len_monomer_subseq < max_hor_len:
                subseq_str = "_".join(subseq)

                if "NM" in subseq_str or "f" in subseq_str:
                    print(subseq_str)
                    exit(-1)
                if subseq_str not in potential_hors:
                    annotation_new_lst = annotation_str.split(subseq_str)
                    annotation_new_str = annotation_str.replace(subseq_str, "X")
                    annotation_new_str = re.sub(r'(X_)\1+', r'\1', (annotation_new_str + "_"))[:-1]

                    new_set_size = len(annotation_new_str.split("_"))
                    potential_hors[subseq_str] = {"set_size": new_set_size, "cnt": len(annotation_new_lst) - 1}
                    potential_hors_names.append(subseq_str)
    for h in potential_hors_all:
        if h not in potential_hors:
            potential_hors_all[h]["set_size"] += len(annotation)
        else:
            potential_hors_all[h]["set_size"] += potential_hors[h]["set_size"]
            potential_hors_all[h]["cnt"] += potential_hors[h]["cnt"]

    for h in potential_hors_names:
        if h not in potential_hors_all:
            potential_hors_names_all.append(h)
            potential_hors_all[h] = {}
            potential_hors_all[h]["set_size"] = potential_hors[h]["set_size"] + set_size
            potential_hors_all[h]["cnt"] = potential_hors[h]["cnt"]
    return potential_hors_all, potential_hors_names_all


def run_iterative_hor_extraction(annotation, known_hors, min_cnt, min_weight, min_hor_len, max_hor_len):
    hors = {}
    hors_lst = []
    h_cnt = 0
    reads = []
    for r in annotation:
        annotation[r] = collapse_annotation(annotation[r])
        reads.append(r)
    reads = sorted(reads)

    hors_log = []
    if len(known_hors) > 0:
        annotation, hors, hors_lst, h_cnt, hors_log = known_hors_annotation(annotation, known_hors, hors, hors_lst, h_cnt, min_cnt)

    while True:
        potential_hors = {}
        potential_hors_names = []
        set_size = 0
        for r in reads:
            annotation_seq = []
            for a in annotation[r]:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            # if h_cnt == 0:
            #     print("_".join(annotation_seq).replace("[1]", ""))
            potential_hors, potential_hors_names = find_potential_hors(annotation[r], annotation_seq, min_hor_len, max_hor_len, hors, potential_hors, potential_hors_names, set_size)
            set_size += len(annotation[r])

        potential_hors_lst = []
        for h in potential_hors_names:
            if potential_hors[h]["cnt"] > min_cnt:
                potential_hors_lst.append([h, potential_hors[h]])
        potential_hors_lst = sorted(potential_hors_lst, key=lambda x: x[1]["set_size"])
        if len(potential_hors_lst) == 0 or set_size - potential_hors_lst[0][1]["set_size"] < min_weight:
            break
        h_cnt += 1
        name = "h" + str(h_cnt)
        hors= build_full_hor(potential_hors_lst[0][0], hors, name)
        hors_lst.append([name, potential_hors_lst[0][0]])
        hors_log.append([name, potential_hors_lst[0][0].replace("[1]",""), hors[name].replace("[1]",""), str(len(hors[name].split("_"))), str(potential_hors_lst[0][1]["cnt"]), \
                                                     str(set_size - potential_hors_lst[0][1]["set_size"]), str(potential_hors_lst[0][1]["set_size"]) ])
        print("\t".join(hors_log[-1]), flush=True)
        for r in reads:
            annotation_seq = []
            for a in annotation[r]:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            annotation[r] = update_annotation(annotation[r], annotation_seq, potential_hors_lst[0][0], name)
            annotation[r] = collapse_annotation(annotation[r])
    return annotation, hors_lst, hors_log

def build_graph(known_hors):
    graph = {}
    graph_rev = {}
    for kh in known_hors:
        chain = kh.replace("[1]", "").split("_")
        for i in range(len(chain)-1):
            if chain[i] not in graph:
                graph[chain[i]] = {}
            graph[chain[i]][chain[i+1]] = 1
            if chain[i+1]+"'" not in graph_rev:
                graph_rev[chain[i+1] + "'"] = {}
            graph_rev[chain[i+1] + "'"][chain[i]+"'"] = 1
    return graph, graph_rev

def run_naive_hor_annotation(annotation, known_hors):
    hors_lst, hors_log = [], []
    kh_graph, kh_graph_rev = build_graph(known_hors)
    new_annotation = {}
    for r in annotation:
        new_annotation[r] = []
        prev_m, prev_cnt, prev_coords = None, None, {}
        for i in range(len(annotation[r])):
            cur_m, cur_cnt, coords = annotation[r][i]
            if (prev_m not in kh_graph and prev_m not in kh_graph_rev) \
                or (prev_m in kh_graph and cur_m not in kh_graph[prev_m]) \
                or (prev_m in kh_graph_rev and cur_m not in kh_graph_rev[prev_m]):
                new_annotation[r].append(annotation[r][i])
            elif (prev_m in kh_graph and cur_m in kh_graph[prev_m]) or (prev_m in kh_graph_rev and cur_m in kh_graph_rev[prev_m]):
                rec = new_annotation[r][-1]
                new_annotation[r][-1] = [rec[0] + "_" + cur_m, 1, {"s": [rec[2]["s"][0]], "e": [coords["e"][-1]] }]
            prev_m, prev_cnt, prev_coords = cur_m, cur_cnt, coords

    hors_mp = {}
    cnt = 0
    for r in new_annotation:
        for i in range(len(new_annotation[r])):
            if "_" in new_annotation[r][i][0]:
                if new_annotation[r][i][0] not in hors_mp:
                    print("h" + str(cnt + 1), new_annotation[r][i][0])
                    hors_mp[new_annotation[r][i][0]] = "h" + str(cnt + 1)
                    cnt += 1
                new_annotation[r][i][0] = hors_mp[new_annotation[r][i][0]]
        new_annotation[r] = collapse_annotation(new_annotation[r])
    return new_annotation, hors_lst, hors_log

def form_hor_dec(annotation, seq):
    new_seq = {}
    for r in reads:
        r_ann, r_seq = annotation[r], seq[r]
        new_seq[r] = []
        i = 0
        for h in r_ann:
            if h[0].startswith("m") or h[0] == "NM":
                for p in range(h[1]):
                    new_seq[r].append(r_seq[i])
                    i += 1
            elif h[0].startswith("h"):
                for p in range(h[1]):
                    sum_idnt = 0
                    mono_h_lst = []
                    j = i
                    while j < len(r_seq) and r_seq[j]["e"] <= h[2]["e"][p]:
                        sum_idnt += r_seq[j]["idnt"]
                        mono_h_lst.append(r_seq[j]["qid"])
                        j += 1
                    start, end = h[2]["s"][p], h[2]["e"][p] #r_seq[i]["s"], r_seq[i + len(mono_h_lst) - 1]["e"]
                    new_seq[r].append({"qid": h[0], "len": len(mono_h_lst), "monomers_lst": mono_h_lst, "s": start, "e": end, "idnt": sum_idnt/len(mono_h_lst)})
                    i += len(mono_h_lst)
    return new_seq

def convert_to_initial_mono(qid, monomers_mp):
    return monomers_mp[qid].split("_")[0]

def convert_to_list_of_monomers(monomers_lst, monomers_mp):
    res = []
    for m in monomers_lst:
        m_name = convert_to_initial_mono(m, monomers_mp)
        res.append(m_name)
    return ",".join(res)

def print_hor_dec(filename, seq, monomers_mp, known_hors_initial):
    prev, prev_qid, start = 0, "", 0
    with open(filename, "w") as fout:
        for r in seq:
            for c in seq[r]:
                if c["qid"] == "NM":
                    if prev_qid != "NM":
                        start = c["s"]
                    prev_qid = "NM"
                else:
                    if prev_qid == "NM":
                        fout.write("\t".join([r, "NM", str(-1), "{0:.2f}".format(55.0), str(start), str(prev), str(prev - start + 1), str(-1)]) + "\n")   
                    if c["qid"] in monomers_mp:
                        fout.write("\t".join([r, convert_to_initial_mono(c["qid"], monomers_mp), str(c["len"]) if not c["qid"].startswith("f") else "-1", "{0:.2f}".format(c["idnt"]), str(c["s"]), str(c["e"]), str(c["e"] - c["s"] + 1), str(c["s"] - prev)]) + "\n")
                    else:
                        fout.write("\t".join([r, convert_to_list_of_monomers(c["monomers_lst"], monomers_mp), str(c["len"]), "{0:.2f}".format(c["idnt"]), str(c["s"]), str(c["e"]), str(c["e"] - c["s"] + 1),  str(c["s"] - prev)]) + "\n")
                prev = c["e"]
                prev_qid = c["qid"]
        if prev_qid == "NM":
            fout.write("\t".join([r, "NM", str(-1), "{0:.2f}".format(55.0), str(start), str(prev), str(prev - start + 1), str(-1)]) + "\n")
    print("Saved to ", filename)


def build_hor_annotation(reads_dec, min_cnt, min_weight, min_len, max_len, monomers_mp, monomers_mp_r, output, canonical, naive_dec):
    annotation = {}
    seq = {}
    for read in reads_dec:
        dec = reads_dec[read]
        seq[read] = []
        annotation[read] = []
        for i in range(len(dec)):
            annotation[read].append([dec[i]["qid"], 1, {"s": [dec[i]["s"]], "e": [dec[i]["e"]], "idnt": dec[i]["idnt"]}])
            seq[read].append({"qid": dec[i]["qid"], "len": 1, "s": dec[i]["s"], "e": dec[i]["e"], "idnt": dec[i]["idnt"]})

    known_hors = []
    known_hors_initial = set()
    if canonical != None:
        print("Canonical HORs identified")
        known_hors, known_hors_initial = build_known_hors(canonical, monomers_mp)

    if naive_dec:
        annotation, hors_lst, hors_log = run_naive_hor_annotation(annotation, known_hors)
    else:
        annotation, hors_lst, hors_log = run_iterative_hor_extraction(annotation, known_hors, min_cnt, min_weight, min_len, max_len)

    for r in annotation:
        annotation_seq = []
        for a in annotation[r]:
            if a[1] > 1:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            else:
                annotation_seq.append(a[0])
        print(r)
        print("_".join(annotation_seq))

    seq = form_hor_dec(annotation, seq)
    print_hor_dec(output, seq, monomers_mp_r, known_hors_initial)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Build HOR decomposition')
    parser.add_argument('sequences', help='fasta-file with long reads or genomic sequences')
    parser.add_argument('monomers', help='fasta-file with monomers')
    parser.add_argument('decomposition', help='tsv-file with monomer decomposition')
    parser.add_argument('output', help='tsv-file to save HOR decomposition')
    parser.add_argument('--canonical', help='txt-file with list of canonical HORs', required = False)
    parser.add_argument('--naive',  help='run naive decomposition using canonical HORs (divides into canonical HORs and their subsequences, --canonical is required)', action="store_true")
    parser.add_argument('--min-idnt',  help='minimum identity of monomer (75 by default)', type=int, default=75, required = False)
    parser.add_argument('--min-reliable',  help='minimum identity of reliable monomer (95, by default)', type=int, default=95, required = False)
    parser.add_argument('--min-cnt',  help='minimum number of potential HOR occurrences to be considered (5 by default)', type=int, default=5, required = False)
    parser.add_argument('--min-weight',  help='minimum weight of potential HOR to be considered (5 by default)', type=int, default=5, required = False)
    parser.add_argument('--min-len',  help='minimum length of HOR in monomers (2 by default)', type=int, default=2, required = False)
    parser.add_argument('--max-len',  help='maximum length of HOR in monomers (30 by default)', type=int, default=30, required = False)
    args = parser.parse_args()

    if args.naive and args.canonical == None:
        print("Naive decomposition requires --canonical to be set")
        exit(-1)

    reads = read_bio_seqs(args.sequences)
    monomers = read_bio_seqs(args.monomers)
    dec, monomers_mp, monomers_mp_r = convert_dec_to_internal_monomers(load_dec(args.decomposition, args.min_idnt, args.min_reliable), monomers)
    filename = args.output

    build_hor_annotation(dec, args.min_cnt, args.min_weight, args.min_len, args.max_len, monomers_mp, monomers_mp_r, args.output, args.canonical, args.naive)


