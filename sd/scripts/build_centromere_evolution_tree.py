#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import numpy as np; np.random.seed(0)
import pandas as pd

import subprocess
import edlib

import graphviz
from graphviz import Digraph

INF = 1000000000
identity_th = 90
cen = "cenX"
tp = "monomer"
freq_id, freq_num = 2, 50
path = "/Sid/tdvorkina/monomers/sdplus_paper/variants/"
cen_fasta = "cenX_0727.fasta"
cen_monomers = "cenX_monomers_hybrids.fasta"
cen_dec = "cenX_decomposition.tsv"
cen_hordec = "cenX_hordecomposition.tsv"
out_dir = "./evolution_tree_result_monomers"

import os
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

def cnt_edist(lst):
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path")
    niceAlign = edlib.getNiceAlignment(result, str(lst[0]), str(lst[1]))
    return niceAlign, result["cigar"], result["editDistance"]

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

def load_dec(filename, reads):
    reads_mapping = {}
    num = 0
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
            if idnt >= identity_th and "H1L" in qseqid:
                reads_mapping[sseqid].append({"sseqid":sseqid, "hor":qseqid.split(".")[0], "qseqid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt, "seq":reads[sseqid].seq[s:e+1], "id":num})
                num += 1
            else:
                reads_mapping[sseqid].append({"sseqid":sseqid, "hor":"?", "qseqid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt, "seq": "", "id": -1})
    return reads_mapping

def load_hordec(filename, reads):
    reads_mapping = {}
    num = 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, mononum, idnt, sstart, send  = ln.strip().split("\t")[:6]
            sseqid = sseqid.split()[0]
            if sseqid not in reads_mapping:
                reads_mapping[sseqid] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            rev = False
            if qseqid.endswith("'"):
                rev = True
                qseqid = qseqid[:-1]
            qseqid = qseqid.split()[0]
            if idnt >= identity_th and "H1L" in qseqid:
                reads_mapping[sseqid].append({"sseqid":sseqid, "qseqid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt, "seq":reads[sseqid].seq[s : e + 1], "id": num, "divseq": [], "mono": []})
                num += 1
    return reads_mapping

def distribute_bymonomers(dec):
    monomers_map = {}
    cnt = 0
    for r in dec:
        for m in dec[r]:
            if m["hor"] != "?":
                if m["qseqid"] not in monomers_map:
                    monomers_map[m["qseqid"]] = []
                monomers_map[m["qseqid"]].append(m)
                cnt += 1
    return monomers_map

def distribute_byhors(hordec, dec):
    for r in hordec:
        i = 0
        for j in range(len(hordec[r])):
            h = hordec[r][j]
            while i < len(dec[r]) and dec[r][i]["e"] < h["s"]:
                i += 1
            while i < len(dec[r]) and dec[r][i]["s"] < h["e"]:
                m = dec[r][i]["qseqid"]
                hordec[r][j]["mono"].append([m, i, dec[r][i]["rev"]])
                i += 1
    return hordec

def distribute_clusters_byhors(hordec, dec):
    i = 0
    aln_id = 0
    for r in hordec:
        for j in range(len(hordec[r])):
            h = hordec[r][j]
            while i < len(dec[r]) and dec[r][i]["e"] < h["s"]:
                i += 1
            while i < len(dec[r]) and dec[r][i]["s"] < h["e"]:
                m = dec[r][i]["qseqid"]
                hordec[r][j]["mono"].append([m, dec[r][i]["id"], dec[r][i]["rev"], dec[r][i]["cluster"], dec[r][i]["s"], dec[r][i]["e"], dec[r][i]["cluster_level"]])
                i += 1
            if hordec[r][j]["qseqid"].endswith("H1L"):
                print(hordec[r][j]["qseqid"], aln_id)
                aln_id += 1
                for m in hordec[r][j]["mono"]:
                    print(" ", m[0], m[3], m[6], m[4]-hordec[r][j]["s"], "-", m[5] - hordec[r][j]["s"], m[1])
    return hordec

def distribute_clusters(dec, clusters):
    i = 0
    for r in dec:
        for j in range(len(dec[r])):
            if i < len(clusters):
                print(j, dec[r][j], i, clusters[i])
            if "id" in dec[r][j] and i < len(clusters) and dec[r][j]["id"] == clusters[i][1]:
                dec[r][j]["cluster"] = clusters[i][2]
                dec[r][j]["cluster_level"] = clusters[i][3]
                i += 1
            else:
                dec[r][j]["cluster"] = []
                dec[r][j]["cluster_level"] = []
    return dec

def extract_horsvs(hordec):
    hor_svs = {}
    for r in hordec:
        for h in hordec[r]:
            if h["qseqid"] not in hor_svs:
                hor_svs[h["qseqid"]] = []
            hor_svs[h["qseqid"]].append(h)
    return hor_svs

def print_isolates(item_id, isolates, num_isolates):
    print("Isolates number:", num_isolates)
    print(" Position, #isolates")
    for i in range(len(isolates)):
        if len(isolates[i]) > 0:
            print(" ", i, len(isolates[i]), isolates[i])

    print("\nSaving isolates distribution to: ", os.path.join(out_dir, cen + "_" + item_id+ "_isolates_bars.png"))
    plt.figure()
    plt.bar([i for i in range(len(isolates))], [len(isolates[i]) for i in range(len(isolates))], color='blue', edgecolor='blue')
    pd_s = pd.Series([len(isolates[i]) for i in range(len(isolates))])
    plt.plot([i for i in range(len(isolates))], pd_s.rolling(window=50, center = True, min_periods = 1).mean().values.tolist(), color='red')
    plt.savefig(os.path.join(out_dir, cen + "_" + item_id + "_isolates_bars.png"))
    plt.close()

def print_runs(all_runs):
    print("\nNumber of runs for f2-f" + str(freq_id) + " < " + str(freq_num) + ": ", len(all_runs))
    print("\nAll runs for f2-f" + str(freq_id) + " < " + str(freq_num) + ": ")
    ii = 0
    for r in sorted(all_runs, key=lambda x: sorted(list(x[0]))[0]):
        print(" ", ii, sorted(list(r[0])), r[1]) 
        ii += 1
    print("\n")
    print("Runs extracted")

def print_tree(n, tree, root, runs, dot, level):
    if n in tree:
        initial_s, initial_e = -1,-1
        if n == root:
            print(n, "root")
            initial_s, initial_e =0, n
        else:
            run = sorted(list(runs[n][0]))
            initial_s, initial_e = run[0], run[-1]
            if len(run) < 13:
                print("id=", n, run, "len=" + str(len(run)), " mutations:", runs[n][1])
            else:
                print("id=", n, "[" + str(run[0]) + "," + str(run[-1]) + "]", "len=" + str(len(run)), " mutations:", runs[n][1])

        runs_sorted = sorted(tree[n], key=lambda x: sorted(list(runs[x][0]))[0])
        for r in runs_sorted:
            run = sorted(list(runs[r][0]))
            s, e, l = run[0], run[-1], len(run)
            if l < 13:
                print(" id=",r, run, "len=" + str(l), " mutations:", runs[r][1])
            else:    
                print(" id=",r, "[" + str(s) + "," + str(e) + "]", "len=" + str(l), " mutations:", runs[r][1])
            dot.node(str(r), str(s) + "-" + str(e))
            dot.edge(str(n), str(r))
        print("")
        for r in runs_sorted:
            print_tree(r, tree, root, runs, dot, level + 1)

def print_arborescence(item_id, tree, runs, num):
    dot = Digraph(comment='Edmonds Result')
    dot.node(str(len(runs)), '0-'+str(num))
    print_tree(len(runs), tree, len(runs), runs, dot, 0)
    dot.render(os.path.join(out_dir, 'test-dot_' + item_id + '_full_INF.gv'))
    if not os.path.exists(os.path.join(out_dir, "subtree_" + item_id + "_dots")):
        os.makedirs(os.path.join(out_dir, "subtree_" + item_id + "_dots"))

    runs_sorted = sorted(tree[len(runs)], key=lambda x: sorted(list(runs[x][0]))[0])
    for r in runs_sorted:
        cur_dot = Digraph(comment='Edmonds Result')
        hors = sorted(list(runs[r][0]))
        s, e, l = hors[0], hors[-1], len(hors)
        cur_dot.node(str(r), str(s) + '-'+str(e))
        if r in tree:
            print_tree(r, tree, len(runs), runs, cur_dot, 0)
            cur_dot.render(os.path.join(out_dir, "subtree_" + item_id + "_dots", "test-dot_" + item_id + "_" + str(s) + "_" + str(e) + "_id" + str(r) + "_INF.gv"))
        else:
            cur_dot.render(os.path.join(out_dir, "subtree_" + item_id + "_dots", "test-dot_" + item_id + "_" + str(s) + "_" + str(e) + "_id" + str(r) + "_empty_INF.gv"))

def print_clusters(item_id, clusters, clusters_depth, clusters_mutations, alns, all_runs):
    cnt = 0
    print("InstanceId, (final) tree nodes it belongs")
    for i in range(len(clusters)):
        print(i, alns[i][0], clusters[i], clusters_mutations[i])
        if len(clusters[i]) > 1:
            cnt += 1
            for it in clusters[i]:
                if len(all_runs[it][0]) < 20:
                    print(" ", sorted(list(all_runs[it][0])), len(all_runs[it][0]), all_runs[it][1])
                else:
                    print(" ", sorted(list(all_runs[it][0]))[0], "-", sorted(list(all_runs[it][0]))[-1], len(all_runs[it][0]), all_runs[it][1])
    print("Total: ", len(clusters), " Potentially Chimeric: ", cnt)

    s, e = 404, 414
    added_nodes = set()
    added_edges = set()
    cur_dot = Digraph(comment='Edmonds Result')
    for j in range(s, e):
        nodes_lst = clusters_depth[j]
        final_nodes_lst = clusters[j]
        print(final_nodes_lst)
        print(nodes_lst)
        prev_n = str(len(all_runs) + j)
        if not prev_n in added_nodes:
            cur_dot.node(str(prev_n), "Monomer=" + str(j), color="red")
            added_nodes.add(prev_n)
        for k in range(len(final_nodes_lst)):
            prev_n = str(len(all_runs) + j)
            n = str(final_nodes_lst[k])
            if not n in added_nodes:
                cur_dot.node(str(n), "node=" + n)
                added_nodes.add(n)
            if not prev_n + "-" + n in added_edges and prev_n != n:
                cur_dot.edge(n, prev_n)
                added_edges.add(prev_n + "-" + n)
            prev_n = n
            for n in nodes_lst[k].split("-"):
                if not n in added_nodes:
                    cur_dot.node(n, "node=" + n)
                    added_nodes.add(n)
                if not prev_n + "-" + n in added_edges and prev_n != n:
                    cur_dot.edge(n, prev_n)
                    added_edges.add(prev_n + "-" + n)
                prev_n = n
    if not os.path.exists(os.path.join(out_dir, "subtree_" + item_id + "_dots")):
        os.makedirs(os.path.join(out_dir, "subtree_" + item_id + "_dots"))
    cur_dot.render(os.path.join(out_dir, "subtree_" + item_id + "_dots",  "test-dot_" + str(s) + "_" + str(e) + "_INF.gv"))


def form_alns(m_instances, m):
    alns = []
    for mi in m_instances:
        if mi["rev"]:
            nice_align, cigar, ed = cnt_edist([mi["seq"], m.seq.reverse_complement()])
        else:
            nice_align, cigar, ed = cnt_edist([mi["seq"], m.seq])
        m_aln = ""
        for i in range(len(nice_align["target_aligned"])):
            if nice_align["target_aligned"][i] != "-":
                m_aln += nice_align["query_aligned"][i]
        alns.append([mi["id"], m_aln])

    for i in range(len(alns)):
        if i > 0 and len(alns[i][1]) != len(alns[i-1][1]):
            print("Alignments don't have equal size: ", alns[i-1][0], alns[i][0])
            print(alns[i-1][1])
            print(alns[i][1])
            exit(-1)
    return alns

def build_freq_map(total_alns):
    monomer_pos = []
    for i in range(len(total_alns[0][1])):
        score = {"A": 0, "C": 0, "G": 0, "T": 0, "-": 0}
        for it in total_alns:
            score[it[1][i]] += 1
        monomer_pos.append(score)
        # pos_scores = sorted([[x, score[x]] for x in score], key=lambda xx: -xx[1])
        # sm = sum([pos_scores[i][1] for i in range(len(pos_scores))])
        # if sm == pos_scores[0][1]:
        #     print(i, pos_scores, "Conserved!")
        # else:
        #     print(i, pos_scores)
    return monomer_pos

# def identify_isolates(item_id, alns, freq_map, distance):
#     isolates = [set() for _ in range(len(alns))]
#     prehistoric_alns = [alns[i][:] for i in range(len(alns))]
#     new_freq_map = []
#     num_isolates = 0
#     for i in range(len(freq_map)):
#         new_freq_map.append(freq_map[i])
#         pos_lst = sorted([[c, freq_map[i][c]] for c in freq_map[i]], key=lambda x: -x[1])
#         for it in pos_lst[1:]:
#             nuc, freq = it
#             left = [-INF for _ in range(len(alns))]
#             nearest = -INF
#             for j in range(len(alns)):
#                 left[j] = nearest
#                 h = alns[j]
#                 if h[1][i] == nuc:
#                     nearest = h[0]
#             right = [INF for _ in range(len(alns))]
#             nearest = INF
#             for j in range(len(alns)-1, -1, -1):
#                 right[j] = nearest
#                 h = alns[j]
#                 if h[1][i] == nuc:
#                     nearest = h[0]
#             for j in range(len(alns)):
#                 h = alns[j]
#                 if h[1][i] == nuc and h[0] - left[j] > distance and right[j] - h[0] > distance:
#                     num_isolates += 1
#                     isolates[j].add(str(i) + "_" + str(nuc))
#                     prehistoric_alns[j][1] = prehistoric_alns[j][1][:i] + pos_lst[0][0] + prehistoric_alns[j][1][i+1:]
#                     new_freq_map[i][pos_lst[0][0]] += 1
#                     new_freq_map[i][nuc] -= 1
#     print_isolates(item_id.replace("/","_"), isolates, num_isolates)
#     return prehistoric_alns, new_freq_map, isolates

def identify_isolates(item_id, alns, freq_map, distance):
    isolates = [set() for _ in range(len(alns))]
    div_pos = [set() for _ in range(len(alns))]
    prehistoric_alns = [alns[i][:] for i in range(len(alns))]
    new_freq_map = []
    num_isolates = 0
    for i in range(len(freq_map)):
        new_freq_map.append(freq_map[i])
        pos_lst = sorted([[c, freq_map[i][c]] for c in freq_map[i]], key=lambda x: -x[1])
        for it in pos_lst[1:]:
            nuc, freq = it
            left = [-100500 for _ in range(len(alns))]
            nearest = -100500
            for j in range(len(alns)):
                left[j] = nearest
                h = alns[j]
                if h[1][i] == nuc:
                    nearest = j
            right = [100500 for _ in range(len(alns))]
            nearest = 100500
            for j in range(len(alns)-1, -1, -1):
                right[j] = nearest
                h = alns[j]
                if h[1][i] == nuc:
                    nearest = j
            for j in range(len(alns)):
                h = alns[j]
                if (h[1][i] == nuc and j - left[j] > distance and right[j] - j > distance):
                    num_isolates += 1
                    isolates[j].add(str(i) + "_" + str(nuc))
                    prehistoric_alns[j][1] = prehistoric_alns[j][1][:i] + pos_lst[0][0] + prehistoric_alns[j][1][i+1:]
                    new_freq_map[i][pos_lst[0][0]] += 1
                    new_freq_map[i][nuc] -= 1

    #print_isolates(item_id.replace("/","_"), isolates, num_isolates)
    return prehistoric_alns, new_freq_map, isolates


def collapse_runs(runs):
    collapsed_runs = []
    used = [False for _ in range(len(runs))]
    for i in range(len(runs)):
        if not used[i]:
            cur_run = [runs[i][0], [runs[i][1][0] ] ]
            used[i] = True
            for j in range(i + 1, len(runs)):
                if len(runs[i][0] & runs[j][0]) == len(runs[i][0]) == len(runs[j][0]):
                    used[j] = True
                    cur_run[1].append(runs[j][1][0])
            collapsed_runs.append(cur_run)
    return collapsed_runs

def extract_runs(alns, freq_map, distance):
    runs = []
    for i in range(len(freq_map)):
        pos_lst = sorted([[c, freq_map[i][c]] for c in freq_map[i]], key=lambda x: -x[1])
        for it in pos_lst[1:freq_id]:
            nuc, freq = it
            if freq < freq_num:
                cur_runs = []
                cur_run = set()
                rightmost = -INF
                for j in range(len(alns)):
                    h = alns[j]
                    if h[1][i] == nuc:
                        if h[0] - rightmost <= distance:
                            cur_run.add(j)
                        else:
                            if len(cur_run) > 0:
                                cur_runs.append(cur_run)
                            cur_run = set()
                            cur_run.add(j)
                        rightmost = h[0]
                if len(cur_run) > 0:
                    cur_runs.append(cur_run)
                for r in cur_runs:
                    runs.append([r, [[nuc, i, freq]] ])
    runs = collapse_runs(runs)
    return runs

def construct_arborescence(item_id, runs, num):
    edges_num = 0
    vertex_num = len(runs) + 1
    gmap = {len(runs):{}}
    for i in range(len(runs)):
        gmap[len(runs)][i] = num - len(runs[i][0])
        edges_num += 1
        for j in range(len(runs)):
            if i != j:
                if runs[j][0].issubset(runs[i][0]):
                    if i not in gmap:
                        gmap[i] = {}
                    gmap[i][j] = len(runs[i][0] - runs[j][0])
                    edges_num +=1
    with open(os.path.join(out_dir, "for_edmond_" + item_id + ".in.graph"), "w") as fout:
        fout.write(str(vertex_num) + " " + str(edges_num) + " " + str(len(runs)) + "\n")
        for n in gmap:
            for n2 in gmap[n]:
                fout.write(str(n+1) + " " + str(n2+1) + " " + str(gmap[n][n2]) + "\n")
    with open(os.path.join(out_dir, "for_edmond_" + item_id + ".out.graph"), 'w') as f:
        #https://github.com/prokls/edmonds-branching-algorithm
        subprocess.run(["python3", "./edmonds.py", os.path.join(out_dir, "for_edmond_" + item_id + ".in.graph")], stdout = f, check = True)
    tree = {}
    with open(os.path.join(out_dir, "for_edmond_" + item_id + ".out.graph"), "r") as fin:
        first_live_line = True
        for ln in fin.readlines():
            if not ln.startswith("c") and not ln.startswith("b"):
                if first_live_line:
                    first_live_line = False
                else:
                    n1, n2 = ln.strip().split()[:2]
                    if int(n1)-1 not in tree:
                        tree[int(n1)-1] = []
                    tree[int(n1)-1].append(int(n2)-1)
    #print_arborescence(item_id, tree, runs, num)
    return tree

def dfs(r, level, level_str, mutations_lst, tree, runs, clusters, clusters_depth, clusters_mutations, prev, min_run_len, maxlevel):
    run = sorted(list(runs[r][0]))
    for h in run:
        if clusters[h][-1] == prev:
            clusters[h][-1] = r
            clusters_depth[h][-1] = level_str
            clusters_mutations[h][-1] = mutations_lst
        else:
            clusters[h].append(r)
            clusters_depth[h].append(level_str)
            clusters_mutations[h].append(mutations_lst)
    if r in tree and level + 1 <= maxlevel:
        for n in tree[r]:
            if len(runs[n][0]) > min_run_len:
                new_mutations_lst = mutations_lst[:]
                dfs(n, level + 1, str(r) + "-" + level_str, new_mutations_lst + runs[n][1], tree, runs, clusters, clusters_depth, clusters_mutations, r, min_run_len, maxlevel)


def construct_clusters(item_id, tree, alns, all_runs):
    res = []
    clusters = [[len(all_runs)] for _ in range(len(alns))]
    clusters_depth = [[str(len(all_runs))] for _ in range(len(alns))]
    clusters_mutations = [[[]] for _ in range(len(alns))]
    min_run_len = 0
    maxlevel = INF
    for r in tree[len(all_runs)]:
        mutations_lst = all_runs[r][1][:]
        dfs(r, 1, str(len(all_runs)), mutations_lst, tree, all_runs, clusters, clusters_depth, clusters_mutations, len(all_runs), min_run_len, maxlevel)
    #print_clusters(item_id.replace("/","_"), clusters, clusters_depth, clusters_mutations, alns, all_runs)
    for i in range(len(alns)):
        res.append([item_id, alns[i][0], clusters[i], clusters_depth[i], clusters_mutations[i]])
    return res

def union_runs(j_mutations, k_mutations, pos):
    res = []
    for p in range(len(j_mutations)):
        if j_mutations[p][0] < pos:
            res.append(j_mutations[p])
        else:
            break
    for p in range(len(k_mutations)):
        if k_mutations[p][0] >= pos:
            res.append(k_mutations[p])
    return res

def are_equal(a_mutations, b_mutations):
    if len(a_mutations) != len(b_mutations):
        return False
    for i in range(len(a_mutations)):
        if a_mutations[i][0] != b_mutations[i][0]:
            return False
        elif set(a_mutations[i][1]) & set(b_mutations[i][1]) == 0:
            return False
    return True

def check_chimerism(i, clusters, alns, distance):
    l, r = max(0, i - distance), min(len(alns), i + distance + 1)
    i_mutations = order_mutations(clusters[i][-1], clusters[i][2])
    for j in range(l, r):
        j_mutations = order_mutations(clusters[j][-1], clusters[j][2])
        if not are_equal(j_mutations, i_mutations):
            for k in range(l, r):
                k_mutations = order_mutations(clusters[k][-1], clusters[k][2])
                if not are_equal(k_mutations, i_mutations):
                    if j != k and j != i and k != i:
                        for pos in range(1, len(alns[j][1]) - 1):
                            chimera_mutations = union_runs(j_mutations, k_mutations, pos)
                            if are_equal(chimera_mutations, i_mutations):
                                print("Monomer", clusters[i][1], " pos=", pos, " Chimera of: Monomer", clusters[j][1], "mutations: ", printable_mutations(j_mutations),\
                                      " and Monomer", clusters[k][1], "mutations:", printable_mutations(k_mutations))
                                return True
    print("Chimera wasn't identified!")
    return False

def order_mutations(mutations, clusters):
    s = []
    for i in range(len(clusters)):
        r = mutations[i]
        cl = clusters[i]
        for mut in r:
            s.append([mut[1], [cl]])
    s = sorted(s, key = lambda x: x[0])
    collapsed_s = []
    if len(s) > 0:
        collapsed_s = [s[0]]
        i = 1
        while i < len(s):
            if collapsed_s[-1][0] == s[i][0]:
                collapsed_s[-1][1].extend(s[i][1])
            else:
                collapsed_s.append(s[i])
            i += 1
    return collapsed_s

def printable_mutations(s):
    res = []
    for it in s:
        res.append(str(it[0]) + "(" + ",".join([str(itt) for itt in it[1]])  + ")")
    return "[ " + ", ".join(res) + " ]"

def process_item(item_id, item, seq, distance):
    print(item_id, distance)
    alns = form_alns(item, seq)
    freq_map = build_freq_map(alns)
    alns, freq_map, isolates = identify_isolates(item_id, alns, freq_map, distance)
    all_runs = extract_runs(alns, freq_map, distance)
    #print_runs(all_runs)
    tree = construct_arborescence(item_id.replace("/", "_"), all_runs, len(alns))
    clusters = construct_clusters(item_id, tree, alns, all_runs)
    print("")
    dist = 20
    print("Distance=", dist)
    for i in range(len(clusters)):
        if len(clusters[i][2]) > 1:
            print(item_id, " Chimeric MonomerId=", clusters[i][1], " mutations", printable_mutations(order_mutations(clusters[i][-1], clusters[i][2])))
        else:
            print(item_id, " MonomerId=", clusters[i][1], " mutations", printable_mutations(order_mutations(clusters[i][-1], clusters[i][2])))
        if len(clusters[i][2]) > 1:
            is_chimera = check_chimerism(i, clusters, alns, dist)
            print(is_chimera)
            print("")
    print("")
    return clusters, tree

def form_consensus(mono_lst, monomers):
    s = ""
    for m in mono_lst:
        if m[0] == "?":
            s += str(["N" for _ in range(171)])
        elif m[2]:
            s += str(monomers[m[0]].seq.reverse_complement())
        else:
            s += str(monomers[m[0]].seq)
    return s

def process_hors(hor_svs, monomers, distance = 12):
    res = []
    for sv in hor_svs:
        if not sv.endswith("H1L"):
            continue
        hor_consensus = SeqRecord(Seq(form_consensus(hor_svs[sv][0]["mono"], monomers)), id=sv, name=sv, description = "")
        clusters, tree = process_item(sv, hor_svs[sv], hor_consensus, distance)
        res.extend(clusters)
    for i in range(len(res)):
        print(res[i])
    return res

def process_monomers(monomers_map, monomers, distance = 12*12):
    res = []
    for m in monomers_map:
        if "/" in m:
            continue
        clusters, tree = process_item(m, monomers_map[m], monomers[m], distance)
        res.extend(clusters)
    res = sorted(res, key = lambda x: x[1])
    exit(-1)
    for i in range(len(res)):
        print(res[i])
    return res

##chimeric HORs begin

def encode(hor, freq_map):
    vec = np.zeros(len(hor))
    for i in range(len(hor)):
        pos_lst = sorted([[c, freq_map[i][c]] for c in freq_map[i]], key=lambda x: -x[1])
        if hor[i] == pos_lst[1][0]:
            vec[i] = 1
    return vec

def mul(v1, v2):
    #res = v1.size - sum(np.logical_xor(v1 == 1, v2 == 1))
    res = int(sum(v1*v2))
    return res

def find_chimeric(alns, freq_map, number=10):
    aln_vecs = []
    for i in range(len(alns)):
        aln_vecs.append(encode(alns[i][1], freq_map))
    for i in range(len(alns)):
        best_j, best_k, best_gain, best_chimera = -1, -1, -INF, -1
        for j in range(max(0, i - number), min(i + number + 1, len(alns))):
            if i != j:
                for k in range(max(0, i - number), min(i + number + 1, len(alns))):
                    if j != k and k != i and j != i:
                        pair_score = max(mul(aln_vecs[i], aln_vecs[j]), mul(aln_vecs[i], aln_vecs[k]))
                        if (aln_vecs[j] == aln_vecs[k]).all():
                            best_p, best_score = 0, pair_score
                        else:    
                            best_p, best_score = -1, -1
                            for p in range(1, len(aln_vecs[j]) - 1):
                                if aln_vecs[j][p] != aln_vecs[k][p]: 
                                    ch_vec = np.concatenate((aln_vecs[j][:p], aln_vecs[k][p:]))
                                    if sum(aln_vecs[i]) >= sum(ch_vec):
                                        ch_score = mul(aln_vecs[i], ch_vec)
                                        if ch_score > best_score:
                                            best_p, best_score = p, ch_score
                            
                        cur_gain = best_score - pair_score
                        if cur_gain > best_gain:
                            best_j, best_k, best_gain, best_chimera = j, k, cur_gain, best_p
        print("")
        print(i, "|B|=", int(sum(aln_vecs[i])), "Chimeric Gain=", best_gain, "sim((A,A'),B)=", max(mul(aln_vecs[i], aln_vecs[best_j]), mul(aln_vecs[i], aln_vecs[best_k])), \
                 "sim(Chimera, B)=", mul(aln_vecs[i], np.concatenate((aln_vecs[best_j][:best_chimera], aln_vecs[best_k][best_chimera:]))), "Best pos=", best_chimera, "Best A", best_k, "Best A'", best_j, flush=True)
    return chimeras

def find_cur_hor_chimeras(item_id, item, seq, distance):
    print(item_id)
    alns = form_alns(item, seq)
    freq_map = build_freq_map(alns)
    alns, freq_map, isolates = identify_isolates(item_id, alns, freq_map, distance)
    find_chimeric(alns, freq_map)

def find_chimeras(hor_svs, monomers, distance = 12):
    for sv in hor_svs:
        if not sv.endswith("H1L"):
            continue
        hor_consensus = SeqRecord(Seq(form_consensus(hor_svs[sv][0]["mono"], monomers)), id=sv, name=sv, description = "")
        find_cur_hor_chimeras(sv, hor_svs[sv], hor_consensus, distance)

##chimeric HORs end

reads = load_fasta(path + cen_fasta, "map")
monomers = load_fasta(path + cen_monomers, "map")
dec = load_dec(path + cen_dec, reads)
hordec = load_hordec(path + cen_hordec, reads)

if tp == "monomer":
    monomers_map = distribute_bymonomers(dec)
    monomer_clusters = process_monomers(monomers_map, monomers)
    dec = distribute_clusters(dec, monomer_clusters)
    hordec = distribute_clusters_byhors(hordec, dec)
elif tp == "hor":
    hordec = distribute_byhors(hordec, dec)
    hor_svs = extract_horsvs(hordec)
    process_hors(hor_svs, monomers)
else:
    hordec = distribute_byhors(hordec, dec)
    hor_svs = extract_horsvs(hordec)
    find_chimeras(hor_svs, monomers)

