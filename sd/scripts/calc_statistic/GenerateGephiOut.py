#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import edlib

c_names = ['black', 'green', 'yellow', 'red', 'firebrick', 'cyan', 'orange', 'chocolate',
          'purple', 'deepskyblue', 'turquoise', 'slategray', 'cornflowerblue', 'blueviolet',
          'violet', 'purple', 'magenta', 'hotpink', 'royalblue', 'azure', 'khaki', 'peru',
           'beige', 'lightgreen', 'thistle', 'oldlace', 'slateblue', 'sienna', 'lime', 'mistyrose',
          'darkturquoise', 'lightcyan', 'olive', 'ghostwhite', 'hotpink', 'tomato', 'darkorange', 'plum', 'lavender']


def cluster_to_color(color):
    return c_names[color + 1]


def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-seqs")
    parser.add_argument("-clusters")
    parser.add_argument("-o")

    return parser.parse_args()


def main():
    args = parse_args()

    f_in = open(args.clusters, 'r')
    colors = []
    sequences = []
    for line in f_in:
        colors.append(int(line.strip().split()[-1]))
    for record in SeqIO.parse(args.seqs, 'fasta'):
        sequences.append(str(record.seq))
    print(np.unique(colors))
    f_in.close()

    seq2 = []
    colors2 = []
    for i in range(len(sequences)):
        isGood = True
        for j in range(0, i):
            if (edlib.align(sequences[i], sequences[j], mode='HW', task='path')['editDistance'] == 0):
                isGood = False
        if isGood:
            seq2.append(sequences[i])
            colors2.append(colors[i])

    sequences = seq2
    colors = colors2

    dist_matrix = np.zeros((len(sequences), len(sequences)), dtype='int')
    for i in range(len(sequences)):
        for j in range(i, len(sequences)):
            dist_matrix[i, j] = edlib.align(sequences[i], sequences[j], mode='HW', task='path')['editDistance']
    dist_matrix = dist_matrix + dist_matrix.T
    print(dist_matrix)
    print(len(sequences))

    edges = np.zeros((len(sequences), len(sequences)), dtype='int')
    edges_list = []
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            if dist_matrix[i, j] < 5:
                edges[i, j] = 5-dist_matrix[i, j]
                edges_list.append((i, j))
    print(len(edges_list))

    f = open(args.o, 'w')
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write(
        '<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n')
    f.write(
        'xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">\n')
    f.write('<key id="d0" for="node" attr.name="color" attr.type="string">\n')
    f.write('<default>yellow</default>\n')
    f.write('</key>\n')
    f.write('<key id="d1" for="edge" attr.name="weight" attr.type="double"/>\n')
    f.write('<graph id="G" edgedefault="undirected">\n')
    for i in range(len(colors)):
        f.write(f'<node id="n{i}">\n')
        f.write(f'<data key="d0">{cluster_to_color(colors[i])}</data>\n')
        f.write('</node>\n')

    for count, (i, j) in enumerate(edges_list):
        f.write(f'<edge id="e{count}" source="n{i}" target="n{j}">\n')
        f.write(f'<data key="d1">{edges[i, j]}</data>\n')
        f.write(f'</edge>')
    f.write('</graph>\n')
    f.write('</graphml>')
    f.close()


if __name__ == "__main__":
    main()