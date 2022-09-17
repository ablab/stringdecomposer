#!/usr/bin/env python3

import argparse
import os
import pathlib
import subprocess
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import edlib
import pandas as pd
import re

from stringdecomposer.py.standard_logger import get_logger
from stringdecomposer.py.git import get_git_revision_short_hash


CUR_FILE = os.path.abspath(__file__)
CUR_DIR = os.path.dirname(CUR_FILE)
SD_BIN = os.path.join(CUR_DIR, 'build', 'bin', 'dp')
LOGREG_FILE = os.path.join(CUR_DIR,
                           'models',
                           'ont_logreg_model.txt')
with open(LOGREG_FILE) as f:
    LR_MODEL_COEF = list(map(float, f.readline().strip().split()))


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


def add_rc_monomers(monomers):
    res = []
    for m in monomers:
        res.append(m)
        res.append(make_record(m.seq.reverse_complement(), m.name + "'", m.id + "'"))
    return res


def convert_to_homo(seq):
    res = ""
    for c in seq:
        if len(res) == 0 or res[-1] != c:
            res += c
    return res


def classify(reads_mapping):
    df = pd.DataFrame(reads_mapping)
    df["idnt_diff"] = df["score"] - df["second_best_score"]
    X = pd.concat([df["score"], df["idnt_diff"]], axis=1, keys = ["idnt", "idnt_diff"])
    X.insert(0, 'intercept', 1)
    y_pred = (X.dot(LR_MODEL_COEF)) > 0
    for i in range(len(reads_mapping)):
        if y_pred[i] != 1:
            reads_mapping[i]["q"] = "?"
    return reads_mapping


def convert_read(decomposition, read, monomers, light = False):
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

    res = classify(res)
    return res


def print_read(fout, fout_alt, dec, read, monomers, identity_th, light):
    dec = convert_read(dec, read, monomers, light)
    for d in dec:
        if d["score"] >= identity_th:
            fout.write("\t".join([read.name, d["m"], d["start"], d["end"], "{:.2f}".format(d["score"]), \
                                                    d["second_best"], "{:.2f}".format(d["second_best_score"]), \
                                                    d["homo_best"], "{:.2f}".format(d["homo_best_score"]), \
                                                    d["homo_second_best"], "{:.2f}".format(d["homo_second_best_score"]), d["q"]]) + "\n")
            for a in d["alt"]:
                star = "-"
                if a == d["m"]:
                    star = "*"
                fout_alt.write("\t".join([read.name, a, d["start"], d["end"], "{:.2f}".format(d["alt"][a]), star]) + "\n")


def convert_tsv(decomposition, reads, monomers, outfile, identity_th, light):
    with open(outfile[:-len(".tsv")] + "_alt.tsv", "w") as fout_alt:
        with open(outfile, "w") as fout:
            cur_dec = []
            prev_read = None
            for ln in decomposition.split("\n")[:-1]:
                read, monomer, start, end = ln.split("\t")[:4]
                read = read.split()[0]
                monomer = monomer.split()[0]
                if read != prev_read and prev_read != None:
                    print_read(fout, fout_alt, cur_dec, reads[prev_read], monomers, identity_th, light)
                    cur_dec = []
                prev_read = read
                start, end = int(start), int(end)
                cur_dec.append({"m": monomer, "start": start, "end": end})
            if len(cur_dec) > 0:
                print_read(fout, fout_alt, cur_dec, reads[prev_read], monomers, identity_th, light)

def run(sequences, monomers, num_threads, scoring, batch_size, raw_file, ed_thr, overlap, logger):
    ins, dels, mm, match = scoring.split(",")
    if not os.path.isfile(SD_BIN):
        logger.info('The binary of String Decomposer is not available. Did you forget to run `make`? Aborting.')
        sys.exit(1)

    logger.info(' '.join(["Run", SD_BIN, "with parameters", sequences, monomers, str(num_threads), str(batch_size), str(overlap), scoring]))
    with open(raw_file, 'w') as f:
        subprocess.run([SD_BIN, sequences, monomers, num_threads, batch_size, overlap, ins, dels, mm, match, str(ed_thr)], stdout = f, check = True)
    with open(raw_file, 'r') as f:
        raw_decomposition = "".join(f.readlines())
    return raw_decomposition



def main():
    parser = argparse.ArgumentParser(description='Decomposes string into blocks alphabet')
    parser.add_argument('sequences', help='fasta-file with long reads or genomic sequences')
    parser.add_argument('monomers', help='fasta-file with monomers')
    parser.add_argument('-t', '--threads',  help='number of threads (by default 1)', default="1", required=False)
    parser.add_argument('-o', '--out-dir',  help='output directory (by default .)', default=".", required=False)
    parser.add_argument('--out-file',  help='output tsv-file (by default "final_decomposition")', default="final_decomposition", required=False)
    parser.add_argument('-i', '--min-identity',  \
                         help='only monomer alignments with percent identity >= MIN_IDENTITY are printed (by default MIN_IDENTITY=0)', type=int, default=0, required=False)
    parser.add_argument('-s', '--scoring', \
                         help='set scoring scheme for SD in the format "insertion,deletion,mismatch,match" (by default "-1,-1,-1,1")', default="-1,-1,-1,1", required=False)
    parser.add_argument('-b', '--batch-size',  help='set size of the batch in parallelization (by default 5000)', type=str, default="5000", required=False)
    parser.add_argument('--second-best', dest="second_best", help='generate second best monomer and homopolymer scores', action="store_true")
    parser.add_argument('--ed_thr', help='align only monomers with edit distance less then ed_thr for each segment (by default align all monomers)', default=-1,
                        type=int, required=False)
    parser.add_argument('-v', '--overlap',  help='set size of batch overlap (by default 500)', type=str, default="500", required=False)
    args = parser.parse_args()
    pathlib.Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    logfn = os.path.join(args.out_dir, 'stringdecomposer.log')
    logger = get_logger(logfn, logger_name='StringDecomposer')

    logger.info(f'cmd: {sys.argv}')
    # TODO get_git_revision_short_hash is commented out
    # since it does not work when stringdecomposer is run from outside of repo
    # logger.info(f'git hash: {get_git_revision_short_hash()}')

    raw_decomp_fn = os.path.join(args.out_dir, args.out_file + "_raw.tsv")
    raw_decomposition = run(args.sequences, args.monomers, args.threads, args.scoring, args.batch_size, raw_decomp_fn, args.ed_thr, args.overlap, logger)
    logger.info("Saved raw decomposition to " + raw_decomp_fn)

    reads = load_fasta(args.sequences, "map")
    monomers = load_fasta(args.monomers)
    monomers = add_rc_monomers(monomers)
    logger.info("Transforming raw alignments...")

    convert_tsv_fn = os.path.join(args.out_dir, args.out_file + ".tsv")
    convert_tsv(raw_decomposition, reads, monomers, convert_tsv_fn, int(args.min_identity), not args.second_best)
    logger.info("Transformation finished. Results can be found in " + convert_tsv_fn)

    logger.info("Thank you for using StringDecomposer!")


if __name__ == "__main__":
    main()
