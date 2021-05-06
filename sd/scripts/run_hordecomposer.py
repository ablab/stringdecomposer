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
import pandas as pd

def load_monodec(filename):
    dec = []
    monomers = set()
    with open(filename, "r") as fin:
        for ln in fin.readlines():
           ref, mon, start, end, idnt = ln.strip().split("\t")[:5]
           dec.append([ref, mon, start, end, idnt])
           monomers.add(mon)
    return dec, monomers

def add_cyclicshifts(hors, hor_name, hor_seq):
    hor_lst = hor_seq.split(",")
    for i in range(len(hor_lst)):
        hors.append([hor_name + '_' + str(i), ",".join(hor_lst[i:] + hor_lst[:i]), len(hor_lst)])
    for i in range(len(hor_lst)):
        hor = hor_lst[i:] + hor_lst[:i]
        hors.append([hor_name + '_' + str(i) + "'", ",".join([it + "'" for it in hor]), len(hor_lst)])
    hor_lst = hor_lst[::-1]
    for i in range(len(hor_lst)):
        hors.append([hor_name + '_' + str(i) + "-", ",".join(hor_lst[i:] + hor_lst[:i]), len(hor_lst)])
    for i in range(len(hor_lst)):
        hor = hor_lst[i:] + hor_lst[:i]
        hors.append([hor_name + '_' + str(i) + "-'", ",".join([it + "'" for it in hor]), len(hor_lst)])
    return hors, len(hor_lst)

def shift(hor_lst):
    min_ind = 0
    for i in range(len(hor_lst)):
        if hor_lst[min_ind] > hor_lst[i]:
            min_ind = i
    return hor_lst[min_ind:] + hor_lst[:min_ind]

def load_horascycle(filename):
    hors, hor_id = {}, {}
    hor_name= ""
    mono_mp = {}
    cnt, hor_cnt = 0, 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if len(ln.split("\t")) < 2:
                continue
            hor_name, hor_seq = ln.strip().split("\t")[:-1]
            hor_lst = hor_seq.split(",")
            hor_lst = shift(hor_lst)
            print(hor_lst)
            mononum = len(hor_lst)
            hor = {}
            for i in range(len(hor_lst) - 1):
                hor[hor_lst[i]] = hor_lst[i + 1]
                hor[hor_lst[i + 1] + "'"] = hor_lst[i] + "'"
            hor[hor_lst[len(hor_lst) - 1]] = hor_lst[0]
            hor[hor_lst[0] + "'"] = hor_lst[len(hor_lst) - 1] + "'"
            hors[hor_name] = [hor, mononum]
            for i in range(len(hor_lst)):
                cnt += 1
                mono_mp[hor_lst[i]] = cnt
                mono_mp[hor_lst[i]+"'"] = -cnt
            cnt += 1
            hor_id[hor_name] = hor_cnt
            hor_cnt += 1
    hors["Mono"] = ["Mono", 1]
    hor_id["Mono"] = 100500
    return hors, hor_id, mono_mp

def load_hors(filename):
    hors = []
    max_len = 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
           hor_name, hor_seq = ln.strip().split("\t")[:-1]
           hors, ln = add_cyclicshifts(hors, hor_name, hor_seq)
           max_len = max(max_len, ln)
    return hors, max_len

def decompose(monodec, hors, max_len):
    idnt_th = 95
    inf = 100500
    dp = [inf for _ in range(len(monodec))]
    dp_rev = ["" for _ in range(len(monodec))]
    dp[0] = 1
    dp_rev[0] = [monodec[0][1], 1, monodec[0][1]]
    for i in range(1, len(dp)):
        set_ln = {}
        for h in hors:
            h_seq, h_len = h[1], h[2]
            if i + 1 < h_len:
               continue
            if h_len in set_ln:
                mono_seq = set_ln[h_len]
            else:
                mono_seq = ",".join([monodec[j][1] for j in range(i - h_len + 1, i + 1)])
            #print(h_seq, mono_seq, dp[i], i - h_len, dp[i - h_len])
            if h_seq == mono_seq:
               if i - h_len < 0 and dp[i] > 1:
                   dp[i] = 1
                   dp_rev[i] = [mono_seq, h_len, h[0]]
               elif i - h_len >= 0 and dp[i] > dp[i - h_len] + 1:
                   dp[i] = dp[i - h_len] + 1
                   dp_rev[i] = [mono_seq, h_len, h[0]]
    hordec = []
    i = len(monodec) - 1
    while i > -1:
       idnt = 0
       for j in range(i - dp_rev[i][1] + 1, i + 1):
           idnt += float(monodec[j][4])
           #print(monodec[j], monodec[j][4])
       idnt /= (i + 1 - (i - dp_rev[i][1] + 1))
       ref, hor, start, end = monodec[i][0], dp_rev[i][0], monodec[i - dp_rev[i][1] + 1][2], monodec[i][3]
       hordec.append([ref, hor, start, end, "{:.2f}".format(idnt), dp_rev[i][2] ])
       i = i - dp_rev[i][1]
    hordec = hordec[::-1]
    return hordec

def decompose_new(monodec, hors, rev=False):
    hordec = []
    if len(monodec) == 0:
        return hordec
    
    fp = 0 if not rev else len(monodec) - 1

    inhor, cur_hor = 1, [monodec[fp]]
    cur_hor_name = None
    for h in hors:
        if monodec[fp][1] in hors[h][0]:
            cur_hor_name = h
    
    bgp, edp, stp = (1, len(monodec), 1) if not rev else (len(monodec) - 2, -1, -1)
    
    for i in range(bgp, edp, stp):
        prev, cur = (cur_hor[-1][1], monodec[i][1]) if stp == 1 else (monodec[i][1], cur_hor[0][1])
#        print(monodec[i], int(cur_hor[-1][3]), monodec[i][2])
        mndist = abs(int(cur_hor[-1][3]) - int(monodec[i][2])) if stp == 1 else abs(int(monodec[i][3]) - int(cur_hor[0][2]))
        if float(monodec[i][4]) > 80 and cur_hor_name in hors and prev in hors[cur_hor_name][0] and hors[cur_hor_name][0][prev] == cur and hors[cur_hor_name][1] > inhor and mndist  < 5:
            inhor += 1
            if stp == 1:
                cur_hor.append(monodec[i])
            else:
                cur_hor = [monodec[i]] + cur_hor
        else:
            idnt = sum([float(x[4]) for x in cur_hor])/len(cur_hor)
            hordec.append([cur_hor[0][0], ",".join([x[1] for x in cur_hor]), cur_hor[0][2], cur_hor[-1][3], "{:.2f}".format(idnt), cur_hor_name ])
            inhor = 1
            cur_hor = [monodec[i]]
            cur_hor_name = "Mono"
            for h in hors:
                if monodec[i][1] in hors[h][0]:
                    cur_hor_name = h
    idnt = sum([float(x[4]) for x in cur_hor])/len(cur_hor)
    hordec.append([cur_hor[0][0], ",".join([x[1] for x in cur_hor]), cur_hor[0][2], cur_hor[-1][3], "{:.2f}".format(idnt), cur_hor_name] )
    return hordec


def handle_one_read(monodec, hors):
    #print("monodec", monodec)
    #print("hors", hors)

    hordec = []
    cur_hor_name = None
    for h in hors:
        if monodec[0][1] in hors[h][0]:
            cur_hor_name = h
    
    def get_hord_pos():
        for i in range(1, len(monodec)):
            prev, cur = monodec[i - 1][1], monodec[i][1]
            if float(monodec[i][4]) < 80 or (cur_hor_name not in hors or prev not in hors[cur_hor_name][0] or hors[cur_hor_name][0][prev] != cur) or (abs(int(monodec[i - 1][3]) - int(monodec[i][2])) > 5):
                   return i

        return 0
   
    hord_pos = get_hord_pos()
    #print("hord pos: ", hord_pos)

    def decsuf(sp):
        md = monodec[sp:]
        return decompose_new(md, hors)

    def decpref(sp):
        md = monodec[:sp]
        return decompose_new(md, hors, rev=True)
    
    hordec = decpref(hord_pos) + decsuf(hord_pos)
    if len(hordec) > 0:
        hordec = hordec[:-1]
    if len(hordec) > 0:
        hordec = hordec[1:]
    
    #print("suffix hord dec: ", hordec)

    return hordec

def decompose_new_reads(monodec, hors):
    hordec = []
    inhor, cur_read_seq, cur_read = 1, [monodec[0]], monodec[0][0]

    md = monodec.copy()
    md.append(["END"])
    for i in range(1, len(md)):
        if cur_read != md[i][0]:
            l, r = 0, len(cur_read_seq) - 1
            while l < len(cur_read_seq) and float(cur_read_seq[l][4]) < 90:
                l += 1
            while r >= 0 and float(cur_read_seq[r][4]) < 90:
                r -= 1
            if l > r:
                cur_read_seq, cur_read = [md[i]], md[i][0]
                continue
            cur_read_seq = cur_read_seq[l: r+1]
            
            hordec += handle_one_read(cur_read_seq, hors)
            cur_read_seq, cur_read = [md[i]], md[i][0]
        else:
            cur_read_seq.append(md[i])
    
    return hordec


def collapse_name(mono_seq, mono_mp, hor, hor_id):
    mono_lst = mono_seq.split(",")
    if len(mono_lst) == 1:
       #if mono_lst[0] not in mono_mp:
       rev = "'" if mono_lst[0].endswith("'") else ""
       if mono_lst[0][1] in "0123456789X":
           return mono_lst[0][0] + rev
       else:
           return mono_lst[0][:2] + rev
       #else:
       #    return str(mono_mp[mono_lst[0]])
    res = ""
    start, end = mono_mp[mono_lst[0]], mono_mp[mono_lst[0]]
    for i in range(1, len(mono_lst)):
        if mono_mp[mono_lst[i-1]] + 1 == mono_mp[mono_lst[i]]:
           end = mono_mp[mono_lst[i]]
        else:
           if end >= start or (start < 0 and end <= start):
              res += str(start) + "-" + str(end) + ","
           else:
              res += str(start) + ","
           start, end = mono_mp[mono_lst[i]], mono_mp[mono_lst[i]]
    if end >= start or (start < 0 and end <= start):
        res += str(start) + "-" + str(end)
    else:
        res += str(start)
    if len(mono_lst) == hor[1] and hor_id == 0:
        if start < 0:
            num = "-" + res[1:].split("-")[0]
        else:
            num = res.split("-")[0]
        return "c<sub>" + num + "</sub>"
    elif len(mono_lst) == hor[1] and hor_id > 0:
        if start < 0:
            num = "-" + res[1:].split("-")[0]
        else:
            num = res.split("-")[0]
        return "c" + str(hor_id) + "<sub>" + num + "</sub>"
    else:
        s, e = res.split("-")[0], res.split("-")[-1]
        if start < 0:
            s = "-" + res[1:].split("-")[0]
        if end < 0:
            e = "-" + res.split("-")[-1]
        return "p<sub>" + s + "-" + e + "</sub>"

def collapse_ivanname(mono_seq, mono_mp, hor):
    mono_lst = mono_seq.split(",")
    if mono_lst[0] not in mono_mp:
        return "hybrid"
    res, rev = mono_mp[mono_lst[0]].split(".")[0] + ".", mono_mp[mono_lst[0]].split(".")[1].endswith("rev")
    if rev:
        cur = int(mono_mp[mono_lst[0]].split(".")[1][:-len("-rev")])
    else:
        cur = int(mono_mp[mono_lst[0]].split(".")[1])
    start, end = cur, cur
    for i in range(1, len(mono_lst)):
        if rev:
           prev, cur = int(mono_mp[mono_lst[i-1]].split(".")[1][:-len("-rev")]), int(mono_mp[mono_lst[i]].split(".")[1][:-len("-rev")])
        else:
           prev, cur = int(mono_mp[mono_lst[i-1]].split(".")[1]), int(mono_mp[mono_lst[i]].split(".")[1])
        if not rev and prev + 1 == cur or rev and prev - 1 == cur:
           end = cur
        else:
           if end != start:
              res += str(start) + "-" + str(end) + "_"
           else:
              res += str(start) + "_"
           start, end = cur, cur
    if end != start:
        res += str(start) + "-" + str(end)
    else:
        res += str(start)
    return res

def collapse_hordec(hordec, mono_mp, hors, hor_ind, ivan_mp):
    hordec_c = [hordec[0]]
    hordec_c_ivan = []
    cnt, idnt = 1, float(hordec[0][4])
    for i in range(1, len(hordec)):
        if hordec[i][1] == hordec_c[-1][1]:
           cnt += 1
           idnt += float(hordec[i][4])
           hordec_c[-1][3] = hordec[i][3]
        else:
           #print(hordec_c[-1][1], collapse_name(hordec_c[-1][1], mono_mp))
           cur_hor_name = hordec_c[-1][5]
           mon_seq = hordec_c[-1][1]
           hordec_c[-1][1] = collapse_name(mon_seq, mono_mp, hors[cur_hor_name], hor_ind[cur_hor_name])
           ivan_name = collapse_ivanname(mon_seq, ivan_mp, hors[cur_hor_name])
           if cnt > 1:
               hordec_c[-1][1] += "<sup>" + str(cnt) + "</sup>"
               ivan_name += "<sup>" + str(cnt) + "</sup>"
           hordec_c[-1][4] = "{:.2f}".format(idnt/cnt)
           hordec_c_ivan.append([x for x in hordec_c[-1]])
           hordec_c_ivan[-1][1] = ivan_name
           idnt, cnt = float(hordec[i][4]), 1
           hordec_c.append(hordec[i])
    hordec_c[-1][4] = "{:.2f}".format(idnt/cnt)
    cur_hor_name = hordec_c[-1][5]
    mon_seq = hordec_c[-1][1]
    hordec_c[-1][1] = collapse_name(mon_seq, mono_mp, hors[cur_hor_name], hor_ind[cur_hor_name])
    ivan_name = collapse_ivanname(mon_seq, ivan_mp, hors[cur_hor_name])
    if cnt > 1:
        hordec_c[-1][1] += "<sup>" + str(cnt) + "</sup>"
        ivan_name += "<sup>" + str(cnt) + "</sup>"
    hordec_c_ivan.append([x for x in hordec_c[-1]])
    hordec_c_ivan[-1][1] = ivan_name
    return hordec_c, hordec_c_ivan

def load_ivan_mapping(filename):
    res = {}
    with open(filename, "r") as fin:
       for ln in fin.readlines():
           mon, ivan_mon = ln.strip().split("\t")
           res[mon] = ivan_mon
           ivan_mon = ivan_mon
           res[mon + "'"] = ivan_mon.split(".")[0] + "'." + ivan_mon.split(".")[1]
    return res


def printStats(hordec_c, outfile):
    cnt_hor = {(oh[-1], 'c' if oh[1][0] == 'c' else oh[1].split('<sup>')[0]): 0 for oh in hordec_c}
    cnt_hor = [[k, sum([1 if '<sup>' not in oh[1] else int(oh[1].split('<sup>')[-1][:-len("</sup>")]) for oh in hordec_c if oh[-1] == k[0] and oh[1].startswith(k[1])])] for k in cnt_hor]
    cnt_hor.sort(key=lambda x: -x[1])
    df = pd.DataFrame(cnt_hor)
    df.to_csv(outfile)


if __name__ == "__main__":
    if len(sys.argv) < 4:
       print("Check arguments! Failed")
       exit(-1)
    monodec, monomers = load_monodec(sys.argv[1])
    #hors, max_len = load_hors(sys.argv[2])
    hors, hor_id, mono_mp = load_horascycle(sys.argv[2])
    outfilename = sys.argv[3]
    ivan_mp = load_ivan_mapping(sys.argv[4])
    #for m in monomers:
    #    hors.append([m, m, 1])
    #hordec = decompose(monodec, hors, max_len)
    hordec = decompose_new_reads(monodec, hors)
    
    with open(outfilename, "w") as fout:
        for i in range(len(hordec)):
           fout.write("\t".join(hordec[i]) + "\n")
    hordec_c, hordec_c_ivan = collapse_hordec(hordec, mono_mp, hors, hor_id, ivan_mp)
    
    printStats(hordec_c, outfilename[:-len(".tsv")] + "stat.csv")

    with open(outfilename[:-len(".tsv")]+"_collapsed.tsv", "w") as fout:
        #for m in sorted(mono_mp.keys()):
        #    fout.write("\t".join([m, str(mono_mp[m])]) + "\n")
        for i in range(len(hordec_c)):
           fout.write("\t".join(hordec_c[i]).replace("{","").replace("}", "") + "\n")
    with open(outfilename[:-len(".tsv")]+"_collapsed_ivannaming.tsv", "w") as fout:
        #for m in sorted(mono_mp.keys()):
        #    fout.write("\t".join([m, str(mono_mp[m])]) + "\n")
        for i in range(len(hordec_c_ivan)):
           fout.write("\t".join(hordec_c_ivan[i]) + "\n")
    complete = 0
    complete_runs = 0
    with open(outfilename[:-len(".tsv")]+"_collapsed.txt", "w") as fout:
#        for m in sorted(mono_mp.keys()):
#            fout.write("\t".join([m, str(mono_mp[m])]) + "\n")
        num = outfilename[len("/Sid/tdvorkina/monomers/HORmon/HORdec/hordec"):].split("_")[0]
        fout.write("<h3>cen" + num + "</h3>\n")
        s = "<p>"
        colors = ["blue", "green", "brown"]
        frequent_partial = {}
        for i in range(len(hordec_c)):
            if hordec_c[i][1].startswith("p"):
               partialhor_name = hordec_c[i][1].split("<sup>")[0]
               if partialhor_name not in frequent_partial:
                   frequent_partial[partialhor_name] = 0
               if "sup" in hordec_c[i][1]:
                   frequent_partial[partialhor_name] += int(hordec_c[i][1].split("<sup>")[1].split("</sup>")[0])
               else:
                   frequent_partial[partialhor_name] += 1
        frequent_partial_lst = sorted([[x, frequent_partial[x]] for x in frequent_partial], key = lambda x: -x[1])
        colored_hors = {}
        print( frequent_partial_lst[:len(colors)])
        for i in range(min(len(colors), len(frequent_partial_lst))):
            colored_hors[frequent_partial_lst[i][0]] = colors[i]
        colored_hors["c"] = "red"
        for i in range(len(hordec_c)):
            hor_name = hordec_c[i][1].split("<sup>")[0]
            if hordec_c[i][1].startswith("c"):
                complete_runs += 1
                s += '<span style="color:red;">' + hordec_c[i][1]+"</span>"
            elif hor_name in colored_hors:
                s += '<span style="color:' + colored_hors[hor_name] + '">' + hordec_c[i][1]+'</span>'
            else:
                s += hordec_c[i][1]
        s += '</p>\n'
        fout.write(s)
    print("Complete hors ", complete)
    print("Complete runs ", complete_runs)
