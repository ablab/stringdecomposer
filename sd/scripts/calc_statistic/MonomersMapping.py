#!/usr/bin/env python

import os
import csv
import numpy as np
import edlib
import math
from subprocess import check_call

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO


path_out = "/Bmo/kolga/runs/SD/cenALLct_v4_IvanA_alphaSat/MonMap/"
path_to_CA_monomers = "/Bmo/kolga/runs/SD/cenALLct_v4_IvanA_alphaSat/cen_mn/"
path_to_Ivan_monomers = "/Bmo/kolga/data/SD/mon_globus/"


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


def getMnPath(cenName):
    return os.path.join(path_to_CA_monomers, cenName, cenName + "_mn.fa")


def getIvanMn(cenName):
    return os.path.join(path_to_Ivan_monomers, cenName + "_mn.fa")


cenName = np.char.array(["cen"]*23) + np.append(np.arange(1, 23).astype(str), "X")
paths_ca = np.char.array(list(map(getMnPath, cenName)))
paths_tsv_res = np.char.array(list(map(lambda x: os.path.join(path_to_CA_monomers, x, "final_decomposition.tsv"), cenName)))
paths_ivans = np.char.array(list(map(getIvanMn, cenName)))
out_dots = np.char.array(list(map(lambda x: os.path.join(path_out, x + "_map.dot"), cenName)))

print(out_dots)

def isH1(x):
    return ("H1" in x.id)

def save_vert(fw, CAmn, IVmn, cenName, cnt_mon):
    fw.write("graph " + cenName + " {\n")
    fw.write(" "* 4 + "rankdir=LR;\n")
    IvFilter = np.array(IVmn)[list(map(isH1, IVmn))]

    for vert in CAmn:
        cnt_mn = 0
        lg = 0.01
        if vert.id in cnt_mon:
            cnt_mn = cnt_mon[vert.id]
            lg = math.log(cnt_mn)
        clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa"]
        print(int(lg))
        fw.write(" " * 4 + '"' + vert.id + "\" [style=filled fillcolor=\"" + clr[int(lg)] + "\" label=\"" + vert.id + "[" + str(cnt_mn) + "]\"];\n")

    for vert in IvFilter:
        fw.write(" " * 4 + '"' + vert.id + "\";\n")


    fw.write(" "*4 + "{rank = same; " + "; ".join(list(map(lambda x: "\"" + x.id + "\"", CAmn))) + ";}\n")
    fw.write(" "*4 + "{rank = same; " + "; ".join(list(map(lambda x: "\"" + x.id + "\"", IvFilter))) + ";}\n")



def get_identity(a, b):
    result = edlib.align(a, b, mode="NW", task="locations")
    if result["editDistance"] == -1:
        return 10 ** 9
    return result["editDistance"] * 100 / max(len(a), len(b))


def rc(seq):
    res_seq = ""
    for i in range(len(seq) - 1, -1, -1):
        if (seq[i] == 'A' or seq[i] == 'a'):
            res_seq += 'T'
        if (seq[i] == 'T' or seq[i] == 't'):
            res_seq += 'A'
        if (seq[i] == 'C' or seq[i] == 'c'):
            res_seq += 'G'
        if (seq[i] == 'G' or seq[i] == 'g'):
            res_seq += 'C'
    return res_seq


def save_edges(fw, CAmn, IVmn, cenN):
    for vt1 in CAmn:
        for vt2 in np.array(IVmn)[np.array(list(map(isH1, IVmn)))]:
            scr = min(get_identity(str(vt1.seq), str(vt2.seq)),get_identity(rc(str(vt1.seq)), str(vt2.seq)))
            if "X" in cenN:
                print(scr)
            thr_wg = [0, 2.5, 5, 11]
            wgs = [5, 3, 1, 0]
            if scr < 10 + 0.001:
                wg = 3
                while (scr < thr_wg[wg]):
                    wg -= 1

                fw.write(" " * 4 + "\"" + vt1.id + "\" -- \""  + vt2.id + "\"")
                fw.write(" [label=\"" + "%.2f" % scr + "\" penwidth="+ str(wgs[wg]) + "];\n")

    fw.write("}\n")

def get_cnt_in_centromers(res_tsv):
    resDiv = 10
    cnt_mon = {}

    with open(res_tsv, "r") as f:
        csv_reader = csv.reader(f, delimiter='\t')
        for row in csv_reader:
            if row[2] == "start":
                continue
            identity = float(row[4])
            if row[-1] == '?':
                continue

            #print("Identity: ", identity)
            if identity < 100 - resDiv:
                continue

            if row[1][-1] == "'":
                row[1] = row[1][:-1]
                continue

            if row[1] not in cnt_mon:
                cnt_mon[row[1]] = 1
            cnt_mon[row[1]] += 1

    return cnt_mon


for i in range(len(cenName)):
    CAmn = load_fasta(paths_ca[i])
    IVmn = load_fasta(paths_ivans[i])
    with open(out_dots[i], "w") as fw:
        cnt_mon = get_cnt_in_centromers(paths_tsv_res[i])
        print(cnt_mon)
        save_vert(fw, CAmn, IVmn, cenName[i], cnt_mon)
        save_edges(fw, CAmn, IVmn, cenName[i])
        #try:
         #   check_call(['dot', '-Tpng', out_dots[i], '-o', out_dots[i][:-3] + "png"])
        #except Exception:
        #    continue