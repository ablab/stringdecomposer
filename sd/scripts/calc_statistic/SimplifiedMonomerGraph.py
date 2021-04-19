#!/bin/usr/env python3
import networkx as nx
from networkx.algorithms import bipartite
from networkx.drawing.nx_agraph import write_dot
from subprocess import check_call
import math
import os
import pandas as pd
from Bio import SeqIO
from SDutils import rc
import SDutils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def GetMaxMatching(kcnt):
    mn_set = {tuple(list(x)[:-1]) for x in kcnt.keys()} | {tuple(list(x)[1:]) for x in kcnt.keys()}

    cnt_mon = {x: 0 for x in mn_set}
    for x in kcnt.keys():
        cnt_mon[tuple(list(x)[:-1])] += kcnt[x]

    B = nx.Graph()
    B.add_nodes_from(list(mn_set), bipartite=0)
    B.add_nodes_from(list([(*x,"_") for x in mn_set]), bipartite=1)

    print("nodes", B.nodes())
    for x, w in kcnt.items():
        B.add_edge(tuple(list(x)[:-1]), tuple((list(x)[1:] + ["_"])), weight=-w)
    for x in mn_set:
        for y in mn_set:
            if list(x)[1:] != list(y)[:-1]:
                B.add_edge(x, (*y, "_"), weight=1000)
                continue

            if (*x, y[-1]) not in kcnt:
                B.add_edge(x, (*y, "_"), weight=0)

    matching = bipartite.matching.minimum_weight_full_matching(B, list(mn_set), "weight")
    return matching


def addNode(G, mn, mncnt, IAnm, IsHybrid):
    lg = 0.01
    if mncnt > 0:
        lg = math.log(mncnt)
    clr = ["red", "#cd5700", "orange", "#ffb300", "yellow", "#ceff1d", "#a7fc00", "#00ff00", "#e0ffff", "#f5fffa"]
    curc = "pink" if IsHybrid else clr[int(lg)]
    G.add_node(mn, style=f'filled', fillcolor=f'{curc}', label=f'{mn}({IAnm}) [{str(int(mncnt))}]')

def addEdges(G, mn1, mn2, cnt, thr):
    scr = cnt
    thr_wg = [100000000, 1000, 500, 100, 1]
    wgs = [7, 5, 3, 1, 0]
    wg = 3
    while (scr > thr_wg[wg]):
        wg -= 1

    if scr > thr:
        G.add_edge(mn1, mn2, label=f'{scr}', penwidth=f'{wgs[wg]}')

def save_seqs(blocks, cluster_seqs_path):
    with open(cluster_seqs_path, "w") as fa:
        for i in range(len(blocks)):
            name = "block" + str(i)
            new_record = SeqRecord(Seq(blocks[i]), id=name, name=name, description="")
            SeqIO.write(new_record, fa, "fasta")


def get_consensus_seq(cluster_seqs_path, arg_threads):
    from Bio.Align.Applications import ClustalwCommandline
    from Bio.Align.Applications import ClustalOmegaCommandline
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Align import MultipleSeqAlignment

    aln_file = '.'.join(cluster_seqs_path.split('.')[:-1]) + "_aln.fasta"
    cmd = ClustalOmegaCommandline(infile=cluster_seqs_path, outfile=aln_file, force=True, threads=arg_threads)
    stdout, stderr = cmd()
    align = AlignIO.read(aln_file, "fasta")

    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.gap_consensus(threshold=0, ambiguous='N')
    consensus = str(consensus).replace('-', '')
    return consensus

def get_blocks(trpl, path_seq, tsv_res):
    blocks = []
    seqs_dict = {}
    for record in SeqIO.parse(path_seq, "fasta"):
        seqs_dict[record.id] = str(record.seq).upper()

    df_sd = pd.read_csv(tsv_res, "\t")
    print(df_sd.head())
    for i in range(1, len(df_sd) - 1):
        if df_sd.iloc[i,4] > 60:
            if df_sd.iloc[i, 1].rstrip("'") == trpl[1]:
                if df_sd.iloc[i - 1, 1].rstrip("'") == trpl[0] and df_sd.iloc[i + 1, 1].rstrip("'") == trpl[2]:
                    blocks.append(seqs_dict[df_sd.iloc[i,0]][df_sd.iloc[i,2]:(df_sd.iloc[i,3] + 1)])
                    if df_sd.iloc[i, 1][-1] == "'":
                        blocks[-1] = rc(blocks[-1])
    return blocks

def SplitMonomers(MnToSplit, mnpath,  sdtsv, path_seq, outd, cenid):
    mnlist = SDutils.load_fasta(mnpath)
    for mn in MnToSplit.keys():
        resmns = [mon for mon in mnlist if mon.id != mn]
        ci = 0
        for ctx in MnToSplit[mn]:
            blocks = get_blocks((ctx[0], mn, ctx[1]), sdtsv, path_seq)
            save_seqs(blocks, os.path.join(outd, "blseq" + cenid + ".fa"))
            consensus = get_consensus_seq(os.path.join(outd, "blseq" + cenid + ".fa"), 16)
            name = mn + "." + str(ci)
            ci += 1
            new_record = SeqRecord(Seq(consensus), id=name, name=name, description="")
            resmns.append(new_record)

        mnlist = resmns

    SDutils.savemn(os.path.join(outd, "u" + cenid + "mn.fa"), mnlist)


def PrintSimplifiedGraph(kcnt, mncnt, CAIA, HybridSet, outd, cenid, path_seq, sdtsv, mnpath, edgeThr=100):
    matching = GetMaxMatching(kcnt)
    G = nx.DiGraph()
    for mn in mncnt.keys():
        addNode(G, mn, mncnt[mn], CAIA[mn], False)

    print(matching)
    cycleid = {mn: i for i, mn in enumerate(list(mncnt.keys()))}
    def redrw(mn1, mn2):
        c1 = cycleid[mn1]
        c2 = cycleid[mn2]

        mnclr = min(c1, c2)
        for mn in cycleid:
            if cycleid[mn] == c1 or cycleid[mn] == c2:
               cycleid[mn] = mnclr

    lwe = set()
    for mn1, mn2 in matching.items():
        if len(mn1) > 1:
            continue

        if (mn1[0], mn2[0]) in kcnt:
            redrw(mn1[0], mn2[0])
            addEdges(G, mn1[0], mn2[0], kcnt[(mn1[0], mn2[0])], edgeThr)
        else:
            lwe.add(mn1[0])
            lwe.add(mn2[0])

    print("CycleID", cycleid)
    for mn1, mn2 in kcnt.keys():
        if cycleid[mn1] != cycleid[mn2] or (mn1 in lwe) or (mn2 in lwe):
            addEdges(G, mn1, mn2, kcnt[mn1, mn2], edgeThr)

    sruno = os.path.join(outd, "simpl" + cenid + ".dot")
    srunopng = os.path.join(outd, "simpl" + cenid + ".png")
    write_dot(G, sruno)
    try:
        check_call(['dot', '-Tpng', sruno, '-o', srunopng])
    except Exception:
        return

    if nx.is_eulerian(G):
        ndCnt = {mn: 0 for mn in mncnt}
        elrCirc = list(nx.eulerian_circuit(G))
        for x, y in elrCirc:
            ndCnt[x] += 1

        print("ElerovCirculatin: ", elrCirc)
        print("MAX vcnt in", max(ndCnt.values()))
        if max(ndCnt.values()) > 1:
            spltNode = {mn: [] for mn, cnt in ndCnt.items() if cnt > 1}
            for i in range(len(elrCirc)):
                x = elrCirc[i][0]
                if ndCnt[x] > 1:
                    spltNode[x].append((elrCirc[i - 1][0], elrCirc[i][1]))
            print("SplitNode: ", spltNode)

            SplitMonomers(spltNode, mnpath, sdtsv, path_seq, outd, cenid)
