from Bio.Seq import Seq
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

sd_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
print(sd_path)
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
                reads_mapping[sseqid].append({"qid": qseqid.split("(")[0], "s": s, "e": e, "rev": rev, "idnt": idnt})
            else:
                reads_mapping[sseqid].append({"qid": "NM", "s": s, "e": e, "rev": False, "idnt": idnt})

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
    monomers_lst = sorted([x.split("(")[0] for x in monomers], key=lambda x: ord(x[0])) #sorted([x for x in monomers], key=lambda x: int(x.split(".")[1].split("/")[0]))
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
        monomers_mp[qid + "'"] = new_qid + "'"
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

def get_annotation_seq_sz(annotation_seq):
    prev = annotation_seq[0]
    cnt = 1
    x_cnt = 0
    for i in range(1, len(annotation_seq)):
        if annotation_seq[i] != prev:
            cnt += 1
            prev = annotation_seq[i]
        if annotation_seq[i] =="X":
            x_cnt +=1
    return cnt, x_cnt

def known_hors_annotation(reads_annotation, known_hors, hors, hors_lst, h_cnt, min_cnt, monomers_mp_r):
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
            set_size += get_annotation_seq_sz(annotation_seq)[0]
            annotation_new_lst = annotation_str.split(kh)
            annotation_new_str = annotation_str.replace(kh, "X")
            #annotation_new_str = re.sub(r'(X_)\1+', r'\1', (annotation_new_str + "_"))[:-1]
            cur_set_size, cur_cnt = get_annotation_seq_sz(annotation_new_str.split("_"))
            new_set_size += cur_set_size
            cnt += cur_cnt
        if cnt == 0:
            continue
        h_cnt += 1
        name = "h" + str(h_cnt)
        hors= build_full_hor(kh, hors, name)
        hors_lst.append([name, kh])
        hors_log.append([name, kh.replace("[1]",""), hors[name].replace("[1]",""), str(len(hors[name].split("_"))), str(cnt), \
                                                     str(set_size - new_set_size), str(new_set_size) ])
        print("\t".join([name, rename(kh.replace("[1]",""), monomers_mp_r), str(len(hors[name].split("_"))), str(cnt), \
                                                     str(set_size - new_set_size), str(new_set_size) ]), flush=True)
        for r in reads_annotation:
            annotation = reads_annotation[r]
            annotation_seq = []
            for a in annotation:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            reads_annotation[r] = update_annotation(annotation, annotation_seq, kh, name)

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
            cnt_sz = sum([annotation[j][2]["sz"] for j in range(i, i + k)])
            new_annotation.append([name, 1, {"s": [annotation[i][2]["s"][0]], "e": [annotation[i + k - 1][2]["e"][-1]], "sz": cnt_sz }])
            i += k
        else:
            new_annotation.append(annotation[i])
            i += 1
    return new_annotation

def update_annotation_superhor(annotation, annotation_seq, new_hor, name):
    new_annotation = []
    k = len(new_hor.split("_"))
    i = 0
    while i < len(annotation):
        if i + k < len(annotation):
            cur_seq = "_".join(annotation_seq[i:i+k])
        else:
            cur_seq = ""
        cur_seq_transformed = re.sub(r"\[\d+\]", "[1]", cur_seq)
        if cur_seq_transformed == new_hor:
            new_annotation.append([name, 1, {"s": [annotation[i][2]["s"][0]], "e": [annotation[i + k - 1][2]["e"][-1]], "sz": k} ])
            i += k
        else:
            new_annotation.append(annotation[i])
            i += 1
    return new_annotation

def collapse_annotation(annotation):
    new_annotation = []
    prev, cnt, cnt_sz, start, end = "", 0, 0, [], []
    for i in range(len(annotation)):
        c, c_cnt = annotation[i][:2]
        if c == prev:
            cnt += c_cnt
            cnt_sz += annotation[i][2]["sz"]
            start.extend(annotation[i][2]["s"])
            end.extend(annotation[i][2]["e"])
        else:
            if prev != "":
                new_annotation.append([prev, cnt, {"s": start, "e": end, "sz": cnt_sz}])
            prev = c
            cnt = c_cnt
            cnt_sz = annotation[i][2]["sz"]
            start = annotation[i][2]["s"]
            end = annotation[i][2]["e"]
    i = len(annotation) - 1
    new_annotation.append([prev, cnt, {"s": start, "e": end, "sz": cnt_sz}])
    return new_annotation

def find_potential_hors(annotation, annotation_seq, min_hor_len, max_hor_len, hors, potential_hors_all, potential_hors_names_all, set_size, superhor):
    potential_hors = {}
    potential_hors_names = []
    annotation_str = "_".join(annotation_seq)
    cur_set_size, _ = get_annotation_seq_sz(annotation_seq)
    for i in range(len(annotation)):
        end_ind = i
        len_subseq, subseq = 0, []
        len_monomer_subseq = 0
        has_diff = False
        has_mono = False
        while end_ind < len(annotation) and len_subseq < max_hor_len \
              and annotation[end_ind][0] != "NM" and not annotation[end_ind][0].startswith("f"):
            len_subseq += annotation[end_ind][2]["sz"]
            subseq.append(annotation_seq[end_ind])
            if annotation[end_ind][0].startswith("m") or superhor:
                has_mono = True
            end_ind += 1
            if min_hor_len < len_subseq < max_hor_len or (min_hor_len < len_subseq and superhor):
                if len(subseq) > 1 and subseq[-1] != subseq[-2]:
                    has_diff = True
                if has_diff and has_mono:
                    subseq_str = "_".join(subseq)
                    if "NM" in subseq_str or "f" in subseq_str:
                        print(subseq_str)
                        exit(-1)
                    if subseq_str not in potential_hors:
                        annotation_new_lst = annotation_str.split(subseq_str)
                        annotation_new_str = annotation_str.replace(subseq_str, "X")
                        #annotation_new_str = re.sub(r'(X_)\1+', r'\1', (annotation_new_str + "_"))[:-1]
                        new_set_size, cnt = get_annotation_seq_sz(annotation_new_str.split("_"))
                        potential_hors[subseq_str] = {"set_size": new_set_size, "cnt": cnt}
                        potential_hors_names.append(subseq_str)
    for h in potential_hors_all:
        if h not in potential_hors:
            potential_hors_all[h]["set_size"] += cur_set_size
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

def rename(st, monomers_mp_r):
    res = ""
    for c in st.split("_"):
        if c.startswith("m"):
            res += monomers_mp_r[c][0]
        else:
            res += c
    return res

def run_iterative_hor_extraction(annotation, known_hors, min_cnt, min_weight, min_hor_len, max_hor_len, monomers_mp_r, superhor = False):
    hors = {}
    hors_lst = []
    h_cnt = 0
    reads = []
    for r in annotation:
        reads.append(r)
    reads = sorted(reads)

    hors_log = []
    if len(known_hors) > 0:
        annotation, hors, hors_lst, h_cnt, hors_log = known_hors_annotation(annotation, known_hors, hors, hors_lst, h_cnt, min_cnt, monomers_mp_r)

    while True:
        potential_hors = {}
        potential_hors_names = []
        set_size = 0
        for r in reads:
            annotation_seq = []
            for a in annotation[r]:
                if superhor:
                    annotation_seq.append(a[0] + "[1]")
                else:
                    annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            # if superhor:
            #     print("_".join(annotation_seq).replace("[1]", ""))
            potential_hors, potential_hors_names = find_potential_hors(annotation[r], annotation_seq, min_hor_len, max_hor_len, hors, potential_hors, potential_hors_names, set_size, superhor)
            set_size += len(annotation[r])

        potential_hors_lst = []
        for h in potential_hors_names:
            if potential_hors[h]["cnt"] > min_cnt:
                potential_hors_lst.append([h, potential_hors[h]])
        potential_hors_lst = sorted(potential_hors_lst, key=lambda x: (x[1]["set_size"], len(x[0].split("_")) ))
        if len(potential_hors_lst) == 0 or set_size - potential_hors_lst[0][1]["set_size"] < min_weight:
            break
        h_cnt += 1
        name = "h" + str(h_cnt)
        hors= build_full_hor(potential_hors_lst[0][0], hors, name)
        hors_lst.append([name, potential_hors_lst[0][0]])
        hors_log.append([name, rename(potential_hors_lst[0][0].replace("[1]",""), monomers_mp_r), str(len(hors[name].split("_"))), str(potential_hors_lst[0][1]["cnt"]), \
                                                     str(set_size - potential_hors_lst[0][1]["set_size"]), str(potential_hors_lst[0][1]["set_size"]) ])
        print("\t".join(hors_log[-1]), flush=True)
        for r in reads:
            annotation_seq = []
            for a in annotation[r]:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            if superhor:
                annotation[r] = collapse_annotation(update_annotation_superhor(annotation[r], annotation_seq, potential_hors_lst[0][0], name))
            else:
                annotation[r] = update_annotation(annotation[r], annotation_seq, potential_hors_lst[0][0], name)
    for r in reads:
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
                new_annotation[r][-1] = [rec[0] + "_" + cur_m, 1, {"s": [rec[2]["s"][0]], "e": [coords["e"][-1]], "sz": len(rec[0].split("_")) + 1}]
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
        #new_annotation[r] = collapse_annotation(new_annotation[r])
    return new_annotation, hors_lst, hors_log

def form_hor_dec(annotation, seq, reads_dec):
    new_seq = {}
    for r in reads_dec:
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

def convert_to_initial_mono_olya(qid, monomers_mp):
    return monomers_mp[qid].split("(")[0]

def convert_to_list_of_monomers_olya(monomers_lst, monomers_mp):
    res = []
    for m in monomers_lst:
        m_name = convert_to_initial_mono_olya(m, monomers_mp)
        res.append(m_name)
    return ",".join(res)

def convert_to_ivan_hors(monomers_lst, monomers_mp, known_hors_initial):
    res = []
    for m in monomers_lst:
        m_name = convert_to_initial_mono(m, monomers_mp)
        res.append(m_name)
    isknown, isknown_rev = False, ""

    if ",".join(res) in known_hors_initial:
        isknown = True
    if ",".join(res[::-1]).replace("'", "") in known_hors_initial:
        isknown = True
        isknown_rev = "'"
    name = []
    prev_hor, prev_num, prev_rev = None, -1, False
    for m in res:
        if m.startswith("S"):
            hor, num = m.split(".")
            rev = False
            if num.endswith("'"):
                rev = True
                num = num[:-1]
            if prev_hor != hor:
                name.append(hor)
                if rev:
                    name.append(num + "-" + num + "'")
                else:
                    name.append(num + "-" + num)
            else:
                if "/" not in num and prev_rev==rev==False and prev_num + 1 == int(num):
                    name[-1] = name[-1].split("-")[0] + "-" + num
                elif "/" not in num and prev_rev==rev==True and prev_num - 1 == int(num):
                    name[-1] = name[-1].split("-")[0] + "-" + num + "'"
                elif rev:
                    name.append(num + "-" + num + "'")
                else:
                    name.append(num + "-" + num)
            if "/" in num:
                prev_hor, prev_num, prev_rev = hor, -100, rev
            else:
                prev_hor, prev_num, prev_rev = hor, int(num), rev
        else:
            name.append(m)
            prev_hor, prev_num, prev_rev = None, -1, False
    new_name = []
    for m in name:
        add = ""
        if m.endswith("'"):
            m = m[:-1]
            add = "'"
        if len(m.split("-")) > 1 and m.split("-")[0] == m.split("-")[1]:
            new_name.append(m.split("-")[0] + add)
        else:
            new_name.append(m + add)
    collapsed_name = []
    prev = None
    for m in new_name:
        if len(collapsed_name) == 0 or prev != m:
            collapsed_name.append(m)
        else:
            m_name = collapsed_name[-1].split("[")[0]
            num = 1
            if len(collapsed_name[-1].split("[")) > 1:
                num = int(collapsed_name[-1].split("[")[1].split("]")[0])
            collapsed_name[-1] = m_name + "[" + str(num + 1) + "]"
        prev = m
    if isknown:
        return collapsed_name[0] + isknown_rev
    else:
        return "_".join(collapsed_name)

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
                        fout.write("\t".join([r, "NM", str(start), str(prev), "{0:.2f}".format(55.0), str(-1), str(prev - start + 1), str(-1)]) + "\n")
                    if c["qid"] in monomers_mp:
                        fout.write("\t".join([r, convert_to_initial_mono_olya(c["qid"], monomers_mp), str(c["s"]), str(c["e"]), "{0:.2f}".format(c["idnt"]),\
                                              str(c["len"]) if not c["qid"].startswith("f") else "-1", \
                                              str(c["e"] - c["s"] + 1), str(c["s"] - prev)]) + "\n")
                    else:
                        fout.write("\t".join([r, convert_to_list_of_monomers_olya(c["monomers_lst"], monomers_mp), str(c["s"]), str(c["e"]), "{0:.2f}".format(c["idnt"]), \
                                              str(c["len"]), str(c["e"] - c["s"] + 1),  str(c["s"] - prev)]) + "\n")
                prev = c["e"]
                prev_qid = c["qid"]
        if prev_qid == "NM":
            fout.write("\t".join([r, "NM", str(-1), "{0:.2f}".format(55.0), str(start), str(prev), str(prev - start + 1), str(-1)]) + "\n")
    print("Saved to ", filename)


def build_hor_annotation(reads_dec, min_cnt, min_weight, min_len, max_len, monomers_mp, monomers_mp_r, output, canonical, naive_dec, run_superhor):
    annotation = {}
    seq = {}
    for read in reads_dec:
        dec = reads_dec[read]
        seq[read] = []
        annotation[read] = []
        for i in range(len(dec)):
            annotation[read].append([dec[i]["qid"], 1, {"s": [dec[i]["s"]], "e": [dec[i]["e"]], "idnt": dec[i]["idnt"], "sz": 1}])
            seq[read].append({"qid": dec[i]["qid"], "len": 1, "s": dec[i]["s"], "e": dec[i]["e"], "idnt": dec[i]["idnt"], "sz": 1})

    known_hors = []
    known_hors_initial = set()
    if canonical != None:
        print("Canonical HORs identified")
        known_hors, known_hors_initial = build_known_hors(canonical, monomers_mp)

    if naive_dec:
        annotation, hors_lst, hors_log = run_naive_hor_annotation(annotation, known_hors)
    else:
        annotation, hors_lst, hors_log = run_iterative_hor_extraction(annotation, known_hors, min_cnt, min_weight, min_len, max_len, monomers_mp_r)

    for r in annotation:
        annotation_seq = []
        for i in range(len(annotation[r])):
            a = annotation[r][i]
            name = a[0]
            # if a[0].startswith("h"):
            #     name = chr(ord("a") + int(a[0][1:]) - 1)
            # else:
            #     name = monomers_mp_r[a[0]][0]
            #     if monomers_mp_r[a[0]]=="NM":
            #         name= "_NM_"
            if a[1] > 1:
                annotation_seq.append(name + "[" + str(a[1]) + "]")
            else:
                annotation_seq.append(name)
        print(r)
        print("".join(annotation_seq))
    if run_superhor:
        print("SuperHORs")
        m_num = len(monomers_mp_r) + 1
        for r in annotation:
            annotation_seq = []
            for i in range(len(annotation[r])):
                a = annotation[r][i]
                name = a[0]
                if a[0].startswith("h"):
                    name = chr(ord("a") + int(a[0][1:]) - 1)
                    if name not in monomers_mp:
                        monomers_mp[name] = "m" + str(m_num)
                        monomers_mp_r["m" + str(m_num)] = name
                        print(name, "m" + str(m_num))
                        m_num += 1
                    annotation[r][i][0] = monomers_mp[name]
                else:
                    name = monomers_mp_r[a[0]][0]
                if a[1] > 1:
                    annotation_seq.append(annotation[r][i][0] + "[" + str(a[1]) + "]")
                else:
                    annotation_seq.append(annotation[r][i][0])
            print(r)
            print("_".join(annotation_seq))
        annotation, hors_lst, hors_log = run_iterative_hor_extraction(annotation, [], min_cnt, min_weight, 1, 1000000000000, monomers_mp_r, True)

        for r in annotation:
            annotation_seq = []
            for a in annotation[r]:
                name = a[0]
                if a[0].startswith("h"):
                    name = a[0] #chr(ord("a") + int(a[0][1:]) - 1)
                else:
                    name = monomers_mp_r[a[0]][0]
                    if monomers_mp_r[a[0]]=="NM":
                        name= "_NM_"
                annotation_seq.append(name)
            print(r)
            print("_".join(annotation_seq))

    seq = form_hor_dec(annotation, seq, reads_dec)
    print_hor_dec(output, seq, monomers_mp_r, known_hors_initial)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Build HOR decomposition')
    parser.add_argument('sequences', help='fasta-file with long reads or genomic sequences')
    parser.add_argument('monomers', help='fasta-file with monomers')
    parser.add_argument('decomposition', help='tsv-file with monomer decomposition')
    parser.add_argument('output', help='tsv-file to save HOR decomposition')
    parser.add_argument('--canonical', help='txt-file with list of canonical HORs', required = False)
    parser.add_argument('--naive',  help='run naive decomposition using canonical HORs (divides into canonical HORs and their subsequences, --canonical is required)', action="store_true")
    parser.add_argument('--superhor',  help='run decomposition into superHORs after naive or classic HOR decomposition', action="store_true")
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

    build_hor_annotation(dec, args.min_cnt, args.min_weight, args.min_len, args.max_len, monomers_mp, monomers_mp_r, args.output, args.canonical, args.naive, args.superhor)


