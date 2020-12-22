import argparse
from collections import Counter
from itertools import groupby
import logging
import os
from string import ascii_uppercase as au
from subprocess import check_call
import sys


import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

SCRIPT_FN = os.path.realpath(__file__)
SCRIPT_DIR = os.path.dirname(SCRIPT_FN)
sys.path.insert(0, os.path.join(SCRIPT_DIR, os.pardir, os.pardir))


from sd.hor.hor_extraction_parser import HORExtractionReport, run_hor_extractor
from sd.sd_parser.sd_parser import run_SD
from sd.standard_logger import get_logger
from sd.utils.bio import read_bio_seqs, calc_identity
from sd.utils.git import get_git_revision_short_hash
from sd.utils.os_utils import cat, expandpath, smart_makedirs
from sd.utils.various import running_mean


logger = logging.getLogger("SD.scripts.align_centromeres")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--asm1", help="Asm1", required=True)
    parser.add_argument("--asm2", help="Asm2", required=True)
    parser.add_argument("-m", "--monomers", help="Monomers", required=True)
    parser.add_argument("--canonical",
                        help='txt-file with a list of canonical HORs',
                        required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--st-asm1", type=int, default=0)
    parser.add_argument("--st-asm2", type=int, default=0)
    parser.add_argument("--match-canonical", type=int, default=1)
    parser.add_argument("--match-noncanonical", type=int, default=20)
    parser.add_argument("--mismatch", type=int, default=0)
    parser.add_argument("--indel", type=int, default=0)
    parser.add_argument("--max-dist", type=int, default=-1)
    parser.add_argument("--trim", type=int, default=5)
    params = parser.parse_args()

    params.asm1 = expandpath(params.asm1)
    params.asm2 = expandpath(params.asm2)
    params.monomers = expandpath(params.monomers)
    params.canonical = expandpath(params.canonical)
    if not os.path.isfile(params.asm1):
        logger.error(f'File does not exists --asm1 == {params.asm1}')
        sys.exit(1)
    if not os.path.isfile(params.asm2):
        logger.error(f'File does not exists --asm2 == {params.asm2}')
        sys.exit(1)
    if not os.path.isfile(params.monomers):
        logger.error(f'File does not exists --monomers == {params.monomers}')
        sys.exit(1)
    if not os.path.isfile(params.canonical):
        logger.error(f'File does not exists --canonical == {params.canonical}')
        sys.exit(1)

    params.outdir = expandpath(params.outdir)
    smart_makedirs(params.outdir)
    return params


def tandem_align(s, t, mc, mr, mm, ind, max_dist, trim):
    logger.info("Alignment starts")
    monosequences = [hor_instance.get_mono_indexes()
                     for hor_instance in s.hor_instances]
    monosequences += [hor_instance.get_mono_indexes()
                      for hor_instance in t.hor_instances]
    common_hor, _ = Counter(monosequences).most_common(1)[0]

    n, m = len(s) + 1, len(t) + 1
    score = [[0 for j in range(m)] for i in range(n)]

    for i in range(1, n):
        score[i][0] = score[i-1][0] - ind
    for j in range(1, m):
        score[0][j] = score[0][j-1] - ind

    for i in range(1, n):
        logger.info(f"Aligning {i} / {n}")
        for j in range(1, m):
            sc = max(score[i-1][j] - ind,
                     score[i][j-1] - ind)
            if s[i-1] is not None and s[i-1] != t[j-1]:
                sc = max(sc,
                         score[i-1][j-1] - mm)
            else:
                if s[i-1] == common_hor:
                    bonus = mc
                else:
                    bonus = mr
                s_nucl = s.get_nucl_segment(i-1)[trim:-trim]
                t_nucl = t.get_nucl_segment(j-1)
                ident, alignment = calc_identity(s_nucl, t_nucl,
                                                 k=max_dist, mode='HW')
                sc = max(sc, score[i-1][j-1] + bonus * ident)
            score[i][j] = sc

    a1, a2 = [], []
    i, j = n-1, m-1
    while i and j:
        if score[i][j] == score[i-1][j] - ind:
            a1.append(s[i-1])
            a2.append('-')
            i -= 1
        elif score[i][j] == score[i][j-1] - ind:
            a1.append('-')
            a2.append(t[j-1])
            j -= 1
        else:
            a1.append(s[i-1])
            a2.append(t[j-1])
            i -= 1
            j -= 1
    if i:
        a1 += s[:i]
        a2 += ['-'] * i
    if j:
        a1 += ['-'] * j
        a2 += t[:j]
    a1.reverse()
    a2.reverse()
    logger.info("Alignment finished")
    return a1, a2


def get_hor_ident(s, t, s_al, t_al, trim, outdir=None, window_size=5):
    logger.info("Calculating HOR alignment identities")
    i, j = 0, 0
    idents = []
    eds = []
    for a, b in zip(s_al, t_al):
        if a == '-':
            j += 1
        elif b == '-':
            i += 1
        else:
            s_nucl = s.get_nucl_segment(i)[trim:-trim]
            t_nucl = t.get_nucl_segment(j)
            ident, alignment = calc_identity(s_nucl, t_nucl, mode='HW')
            ed = alignment['editDistance']
            idents.append(ident)
            eds.append(ed)
            i += 1
            j += 1
    assert i == len(s) and j == len(t)
    if outdir is not None:
        rm = running_mean(idents, window_size=window_size)
        outtxtfn = os.path.join(outdir, 'hor_alignment_ident.txt')
        logger.info(f'Saving HOR identities to {outtxtfn}')
        with open(outtxtfn, 'w') as f:
            for i, v in enumerate(rm):
                print(i, v, file=f)
        plt.figure(figsize=(15, 5), dpi=300)
        plt.plot(rm)
        plt.ylim(0.9, 1)
        plt.title('HOR alignment identities', fontsize=20)
        plt.xlabel('HOR', fontsize=15)
        plt.ylabel('Identity', fontsize=15)
        outpdffn = os.path.join(outdir, 'hor_alignment_ident.pdf')
        logger.info(f'Plotting HOR identities to {outpdffn}')
        plt.savefig(outpdffn, format='pdf')
        plt.close()
    return idents, eds


def alignment2substrs(s_al, t_al):
    i, j = 0, 0
    align = []
    for a, b in zip(s_al, t_al):
        if a == '-':
            j += 1
        elif b == '-':
            i += 1
        else:
            align.append((i, j))
            i += 1
            j += 1
    substrs = []
    s1, s2 = align[0]

    s = [a for a in s_al if a != '-']
    t = [b for b in t_al if b != '-']
    for (i0, j0), (i1, j1) in zip(align[0:], align[1:]):
        assert s[i0] == t[j0]
        if i1 - i0 == 1 and j1 - j0 == 1:
            continue
        assert s[s1:i0+1] == t[s2:j0+1]
        substrs.append(((s1, i0+1), (s2, j0+1)))
        s1, s2 = i1, j1
    if i1 != s1:
        assert s[s1:i1+1] == t[s2:j1+1]
        substrs.append(((s1, i1+1), (s2, j1+1)))
    return substrs


def plot_shared_HOR_runs(substrs, str1, str2, outdir, min_len=1):
    fig, ax = plt.subplots(2, figsize=(25, 8), dpi=300)
    substrs = filter(lambda x: x[0][1]-x[0][0] >= min_len, substrs)
    substrs = list(substrs)

    color = plt.cm.rainbow(np.linspace(0, 1, len(substrs)))
    np.random.seed(0)
    np.random.shuffle(color)

    for ((a, b), (c, d)), col in zip(substrs, color):
        rect1 = matplotlib.patches.Rectangle((a, 0),
                                             b-a, 1,
                                             color=col,
                                             alpha=0.7)
        rect2 = matplotlib.patches.Rectangle((c, 0),
                                             d-c, 1,
                                             color=col,
                                             alpha=0.7)

        ax[0].add_patch(rect1)
        ax[1].add_patch(rect2)

    xmax = max(substrs[-1][0][1], substrs[-1][1][1])
    ax[0].set_xlim(0, xmax)
    ax[1].set_xlim(0, xmax)
    ax[1].set_xticks(np.arange(0, xmax, 100))
    ax[1].set_xticklabels(labels=np.arange(0, xmax, 100),
                          rotation=30,
                          fontsize=14)

    ax[0].tick_params(
            axis='x',           # changes apply to the x-axis
            which='both',       # both major and minor ticks are affected
            bottom=False,       # ticks along the bottom edge are off
            top=False,          # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off

    ax[0].tick_params(
            axis='y',
            which='both',
            left=False,
            right=False,
            labelleft=False)
    ax[1].tick_params(
            axis='y',
            which='both',
            left=False,
            right=False,
            labelleft=False)
    ax[0].set_ylabel(str1.seq_id[:10], fontsize=20)
    ax[1].set_ylabel(str2.seq_id[:10], fontsize=20)
    outfn = os.path.join(outdir, 'shared_hor_runs.pdf')
    logger.info(f'Plotting shared HOR runs in {outfn}')
    plt.savefig(outfn, format='pdf')
    plt.close()


def translate_hor_ms(a):
    if a != '-' and a is not None:
        if len(a) > 1:
            trans_a = f'{au[a[0]]}'
            last = 0
            for i in range(1, len(a)):
                if a[i] != a[i-1] + 1:
                    if last == i-1:
                        trans_a += f'_{au[a[i]]}'
                    else:
                        trans_a += f'-{au[a[i-1]]}_{au[a[i]]}'
                    last = i
            trans_a += f'-{au[a[-1]]}'
            a = trans_a
        else:
            a = au[a[0]]
    return a


def compress_alignments(a1, a2, asm1_name, asm2_name, outfn):
    logger.info(f'Writing human readable compressed alignments to {outfn}')
    with open(outfn, 'w') as f:
        print(asm1_name, asm2_name, 'cnt', sep='\t', file=f)
        als = list(zip(a1, a2))
        for p, g in groupby(als):
            x = len(list(g))
            a, b = p
            a = translate_hor_ms(a)
            b = translate_hor_ms(b)
            print(a, b, x, sep='\t', file=f)


def substr2coords(substr, s, t, min_len=10, outdir=None,
                  shift_s=0, shift_t=0):
    coords = []
    for (s_st, s_en), (t_st, t_en) in substr:
        s_len = s_en - s_st
        t_len = t_en - t_st
        if s_len <= min_len or t_len <= min_len:
            continue
        s_n_st = shift_s + s.get_hor_instance_coords(s_st,   original=True)[0]
        s_n_en = shift_s + s.get_hor_instance_coords(s_en-1, original=True)[1]
        t_n_st = shift_t + t.get_hor_instance_coords(t_st,   original=True)[0]
        t_n_en = shift_t + t.get_hor_instance_coords(t_en-1, original=True)[1]
        if s_n_en < s_n_st and t_n_st < t_n_en:
            s_n_st, s_n_en = s_n_en + 1, s_n_st + 1
            t_n_st, t_n_en = t_n_en + 1, t_n_st + 1
        coords.append((s_n_st, s_n_en, t_n_st, t_n_en))
    if outdir is not None:
        outfn = os.path.join(outdir, f'shared_coords_minLen{min_len}.tsv')
        with open(outfn, 'w') as f:
            print('st|' + s.seq_id, 'en|' + s.seq_id,
                  'st|' + t.seq_id, 'en|' + t.seq_id,
                  sep='\t', file=f)
            for s_n_st, s_n_en, t_n_st, t_n_en in coords:
                print(s_n_st, s_n_en, t_n_st, t_n_en, sep='\t', file=f)
    return coords


def alignment_graph(s_al, t_al, outdir=None):
    def compress_alignment(als):
        compr = []
        for p, g in groupby(als):
            x = len(list(g))
            p = translate_hor_ms(p)
            compr.append((p, x))
        compr = ','.join(f'{p}[{x}]' for p, x in compr)
        compr = f' {compr} '
        return compr

    substrs = []
    assert len(s_al) == len(t_al)
    i, j = 0, 0
    while i < len(s_al):
        if s_al[i] == t_al[i]:
            while j < len(s_al) and s_al[j] == t_al[j]:
                j += 1
            mode = 'm'
            substr = compress_alignment(s_al[i:j])
        elif s_al[i] == '-':
            while j < len(s_al) and s_al[j] == '-':
                j += 1
            mode = '-*'
            substr = compress_alignment(t_al[i:j])
        elif t_al[i] == '-':
            while j < len(t_al) and t_al[j] == '-':
                j += 1
            mode = '*-'
            substr = compress_alignment(s_al[i:j])
        else:
            while j < len(t_al) and s_al[j] != t_al[j]:
                j += 1
            substr = (compress_alignment(s_al[i:j]),
                      compress_alignment(t_al[i:j]))
            mode = 'mm'
        substrs.append((i, j, substr, mode))
        i = j

    graph = nx.MultiDiGraph()
    # graph.graph['graph']={'rankdir':'LR', fontsize=12}
    node = 0
    for _, _, substr, mode in substrs:
        if mode == 'm':
            graph.add_edge(node, node + 1, label=substr, color='green')
            node += 1
        elif mode == 'mm':
            graph.add_edge(node, node + 1, label=substr[0], color='blue')
            graph.add_edge(node, node + 1, label=substr[1], color='red')
            node += 1
        elif mode == '-*':
            graph.add_edge(node, node, label=substr, color='red')
        elif mode == '*-':
            graph.add_edge(node, node, label=substr, color='blue')

    for node in graph.nodes:
        graph.nodes[node]['label'] = ''
        graph.nodes[node]['width'] = 0.1
        graph.nodes[node]['height'] = 0.1

    if outdir is not None:
        dotfn = os.path.join(outdir, 'alignment_graph.dot')
        pdffn = os.path.join(outdir, 'alignment_graph.pdf')
        nx.drawing.nx_pydot.write_dot(graph, dotfn)
        cmd = ['dot', '-Tpdf', dotfn, '-o', pdffn]
        check_call(cmd)
    return graph


def align_centromeres(asm1_fn, asm2_fn,
                      monomers_fn,
                      outdir,
                      canonical_fn,
                      n_threads,
                      match_canonical,
                      match_noncanonical,
                      mismatch,
                      indel,
                      max_dist,
                      trim,
                      st_asm1,
                      st_asm2):
    asm1_name = list(read_bio_seqs(asm1_fn).keys())[0]
    asm2_name = list(read_bio_seqs(asm2_fn).keys())[0]
    logger.info(f'Assembly 1 name = {asm1_name}')
    logger.info(f'Assembly 2 name = {asm2_name}')
    sd_asm1_fn = run_SD(sequences_fn=asm1_fn,
                        monomers_fn=monomers_fn,
                        outdir=os.path.join(outdir,
                                            f'SD_report_{asm1_name}'),
                        n_threads=n_threads)
    sd_asm2_fn = run_SD(sequences_fn=asm2_fn,
                        monomers_fn=monomers_fn,
                        outdir=os.path.join(outdir,
                                            f'SD_report_{asm2_name}'),
                        n_threads=n_threads)
    sd_asm_fn = os.path.join(outdir, 'SD_report_concat.tsv')
    logger.info(f'Cat {sd_asm1_fn} and {sd_asm2_fn} into {sd_asm_fn}')
    cat([sd_asm1_fn, sd_asm2_fn], sd_asm_fn)

    asm_fn = os.path.join(outdir, 'asm.fasta')
    logger.info(f'Cat {asm1_fn} and {asm2_fn} into {asm_fn}')
    cat([asm1_fn, asm2_fn], asm_fn)

    hor_report_fn = run_hor_extractor(sequences_fn=asm_fn,
                                      sd_report_fn=sd_asm_fn,
                                      monomers_fn=monomers_fn,
                                      canonical_fn=canonical_fn,
                                      outdir=os.path.join(outdir,
                                                          'hor_dec'))
    hor_report = HORExtractionReport(hor_report_fn=hor_report_fn,
                                     sd_report_fn=sd_asm_fn,
                                     monomers_fn=monomers_fn,
                                     sequences_fn=asm_fn)

    hor_asm1 = hor_report.horstring_set[asm1_name]
    hor_asm2 = hor_report.horstring_set[asm2_name]

    align1, align2 = tandem_align(hor_asm1, hor_asm2,
                                  mc=match_canonical,
                                  mr=match_noncanonical,
                                  mm=mismatch,
                                  ind=indel,
                                  max_dist=max_dist,
                                  trim=trim)
    ident, hor = get_hor_ident(s=hor_asm1, t=hor_asm2,
                               s_al=align1, t_al=align2,
                               trim=trim,
                               outdir=outdir)

    substr = alignment2substrs(align1, align2)
    plot_shared_HOR_runs(substr, hor_asm1, hor_asm2, outdir)
    compress_alignments(align1, align2, asm1_name, asm2_name,
                        outfn=os.path.join(outdir, 'human_readable_align.tsv'))

    for min_len in [1, 5, 10, 20, 50, 100]:
        substr2coords(substr=substr,
                      s=hor_asm1,
                      t=hor_asm2,
                      min_len=min_len,
                      outdir=outdir,
                      shift_s=st_asm1, shift_t=st_asm2)
    alignment_graph(s_al=align1, t_al=align2, outdir=outdir)


def main():
    params = parse_args()
    logfn = os.path.join(params.outdir, 'align_centromeres.log')
    global logger
    logger = get_logger(logfn,
                        logger_name='SD')
    logger.info(f'Aligning Centromeres with {SCRIPT_FN} started')
    logger.info('cmd: {}'.format(sys.argv))
    logger.info('git hash: {}'.format(get_git_revision_short_hash()))

    align_centromeres(asm1_fn=params.asm1,
                      asm2_fn=params.asm2,
                      monomers_fn=params.monomers,
                      outdir=params.outdir,
                      canonical_fn=params.canonical,
                      n_threads=params.threads,
                      match_canonical=params.match_canonical,
                      match_noncanonical=params.match_noncanonical,
                      mismatch=params.mismatch,
                      indel=params.indel,
                      max_dist=params.max_dist,
                      trim=params.trim,
                      st_asm1=params.st_asm1,
                      st_asm2=params.st_asm2)


if __name__ == "__main__":
    main()
