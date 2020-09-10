#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import os
from os import listdir
from os.path import join, isfile
import sys
import argparse

import subprocess

import numpy as np; np.random.seed(0)
import pandas as pd

import re
import edlib

import joblib


def edist(lst):
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
    ed, cigar = edist([str(p1), str(p2)])
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


def classify(reads_mapping, clf):
    df = pd.DataFrame(reads_mapping)
    df["idnt_diff"] = df["score"] - df["second_best_score"]
    X = pd.concat([df["score"], df["idnt_diff"]], axis=1, keys = ["idnt", "idnt_diff"])
    X_scaled = X
    y_pred = list(clf.predict(X_scaled))
    for i in range(len(reads_mapping)):
        if y_pred[i] != 1:
            reads_mapping[i]["q"] = "?"
    return reads_mapping


def convert_to_homo(seq):
    res = ""
    for c in seq:
        if len(res) == 0 or res[-1] != c:
            res += c
    return res


def convert_read(decomposition, read, monomers, clf, light = False):
    res = []
    for d in decomposition:
        monomer, start, end = d["m"], d["start"], d["end"]
        if light:
            scores = {}
            for m in monomers:
                if m.name == monomer:
                    score = aai([read.seq[start:end + 1], m.seq])
                    scores[m.name] = score
            res.append({"m": monomer, "start": str(d["start"]), "end": str(d["end"]), "score": scores[monomer], \
                                    "second_best": "None", "second_best_score": -1,\
                                    "homo_best": "None", "homo_best_score": -1,\
                                    "homo_second_best": "None", "homo_second_best_score": -1,\
                                    "alt": {}, "q": "+"})
        else:
            scores = {}
            for m in monomers:
                score = aai([read.seq[start:end + 1], m.seq])
                scores[m.name] = score
            if monomer == None:
                for s in scores:
                    if monomer == None or scores[s] > scores[monomer]:
                        monomer = s
            secondbest, secondbest_score = None, -1
            for m in scores:
                if m != monomer: # and abs(scores[m] - scores[monomer]) < 5:
                    if not secondbest or secondbest_score < scores[m]:
                        secondbest, secondbest_score = m, scores[m]

            homo_scores = []
            homo_subseq = convert_to_homo(read.seq[start:end + 1])
            for m in monomers:
                score = aai([homo_subseq, convert_to_homo(m.seq)])
                homo_scores.append([m.name, score])
            homo_scores = sorted(homo_scores, key = lambda x: -x[1])
            res.append({"m": monomer, "start": str(d["start"]), "end": str(d["end"]), "score": scores[monomer], \
                                    "second_best": str(secondbest), "second_best_score": secondbest_score,\
                                    "homo_best": homo_scores[0][0], "homo_best_score": homo_scores[0][1],\
                                    "homo_second_best": homo_scores[1][0], "homo_second_best_score": homo_scores[1][1],\
                                    "alt": scores, "q": "+"})

    res = classify(res, clf)
    return res
