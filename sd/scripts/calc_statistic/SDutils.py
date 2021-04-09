import os
import sys
import edlib
from Bio import SeqIO
import pandas as pd
import numpy as np

pathToSD = "/Bmo/kolga/soft/stringdecomposer/sd/run_decomposer.py"

def get_dist(a, b):
    result = edlib.align(a, b, mode="NW", task="locations")
    return result["editDistance"]


def seq_identity(a, b):
    result = edlib.align(a, b, mode="NW", task="locations")
    if result["editDistance"] == -1:
        return 10**9
    return result["editDistance"] * 100 / max(len(a), len(b))


def get_closest_mn(mn, mn_list, distf=get_dist):
    dists = [(distf(str(mn.seq), str(mnx.seq)), mnx.id) for mnx in mn_list]
    return min(dists)


def map_mn(mn_list1, mn_list2, distf=get_dist):
    return {mn1.id: get_closest_mn(mn1, mn_list2, distf) for mn1 in mn_list1}


def get_blocks(path_seq, tsv_B_res, mnidList=None):
    seqs_dict = {}
    for record in SeqIO.parse(path_seq, "fasta"):
        seqs_dict[record.id] = str(record.seq).upper()

    df_sd = pd.read_csv(tsv_B_res, "\t")
    blocks = []
    for i in range(len(df_sd)):
        if df_sd.iloc[i,4] > 60:
            if mnidList is None or df_sd.iloc[i, 1].rstrip("'") in mnidList:
                blocks.append(seqs_dict[df_sd.iloc[i,0]][df_sd.iloc[i,2]:(df_sd.iloc[i,3] + 1)])
                if df_sd.iloc[i, 1][-1] == "'":
                    blocks[-1] = rc(blocks[-1])
    return blocks

def get_monocent(tsv_res):
    df_sd = pd.read_csv(tsv_res, "\t")
    mncent = []
    for i in range(len(df_sd)):
        if df_sd.iloc[i, 4] > 60:
            if "cen1_" in df_sd.iloc[i, 0] and df_sd.iloc[i, 1][-1] == "'":
                continue
            mncent.append(df_sd.iloc[i, 1])
    return mncent


def rc(blc):
    change = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
    return "".join([change[ch] for ch in blc.upper()[::-1]])


def process_readline(line, is_python3=sys.version.startswith("3.")):
    if is_python3:
        return str(line, "utf-8").rstrip()
    return line.rstrip()

def sys_call(cmd):
    import shlex
    import subprocess

    if isinstance(cmd, list):
        cmd_list = cmd
    else:
        cmd_list = shlex.split(cmd)

    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    while not proc.poll():
        line = process_readline(proc.stdout.readline())
        if line:
            print(line)
        if proc.returncode is not None:
            break

    for line in proc.stdout.readlines():
        line = process_readline(line)
        if line:
            print(line)

    if proc.returncode:
        print("system call for: \"%s\" finished abnormally, OS return value: %d" % (cmd, proc.returncode))

def unique(lst_rec):
    return [val for key, val in {x.id : x for x in lst_rec}.items()]

def load_fasta(filename):
    return list(SeqIO.parse(filename, "fasta"))

def run_SD(pathToMon, seqPath, outName):
    os.makedirs(outName, exist_ok=True)
    pathToMon =  os.path.abspath(pathToMon)
    seqPath =  os.path.abspath(seqPath)
    origin_dir = os.path.abspath(os.getcwd())
    os.chdir(outName)
    if not os.path.exists("final_decomposition.tsv"):
        sys_call(["python3", pathToSD, seqPath, pathToMon, "-t", "30", "--fast"])
    os.chdir(origin_dir)
    return os.path.join(outName, "final_decomposition.tsv")


def get_sq_sum(blocks, dists):
    dists = np.array(dists)
    dists = dists.min(axis=1)
    dists **=2
    return (sum(dists)/len(blocks))**0.5


def blocks_sqmean(path_seq, tsv, mns):
    blocks = get_blocks(path_seq, tsv)
    dists = [[min(seq_identity(str(mn.seq), bl), seq_identity(str(mn.seq), rc(bl))) for mn in mns] for bl in blocks]
    sq_sum = get_sq_sum(blocks, dists)
    return sq_sum


def DaviesBouldinIndex(path_seq, tsv, mns):
    blocks = get_blocks(path_seq, tsv)
    dists = [[min(seq_identity(str(mn.seq), bl), seq_identity(str(mn.seq), rc(bl))) for mn in mns] for bl in blocks]
    blmn = [dists[i].index(min(dists[i])) for i in range(len(blocks))]
    radius = [0] * len(mns)
    for i in range(len(blocks)):
        radius[blmn[i]] = max(radius[blmn[i]], dists[i][blmn[i]])

    DBindex = 0
    for i in range(len(mns)):
        dbmx = 0
        for j in range(len(mns)):
            if i == j:
                continue
            if seq_identity(mns[i].seq, mns[j].seq) == 0:
                dbmx = 100
                continue
            if seq_identity(mns[i].seq, rc(mns[j].seq)) == 0:
                dbmx = 100
                continue

            #print(seq_identity(mns[i].seq, mns[j].seq))
            #print(seq_identity(mns[i].seq, rc(mns[j].seq)))
            dbmx = max(dbmx,
                       (radius[i] + radius[j])/min(seq_identity(mns[i].seq, mns[j].seq),
                                                   seq_identity(mns[i].seq, rc(mns[j].seq))))

        DBindex += dbmx
    DBindex /= len(mns)
    return DBindex


def savemn(out, mns):
    with open(out, "w") as fa:
        for record in mns:
            SeqIO.write(record, fa, "fasta")