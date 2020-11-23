#!/usr/bin/env python

import os
import sys
import csv
import argparse
import shutil

def parse_args():
    parser = argparse.ArgumentParser(description='Merge String Decomposer output and Ivans monomers annotation')
    parser.add_argument('iAnnot', help='path to csv fila with Ivans monomers annotaion')
    parser.add_argument('SDout', help='path to  final_decomposition.tsv')
    parser.add_argument('o', help="path to output file")
    return parser.parse_args()

def main():
    cen_start = 57800000
    args = parse_args()

    #Contig, Start, End, Name, Score, Starnd, Length, Gap
    iAnnotList = []
    with open(args.iAnnot, 'r') as ifl:
        csv_reader = csv.reader(ifl)
        for line in csv_reader:
            if line[0] != "Contig":
                iAnnotList.append(line)
                iAnnotList[-1][1] = int(iAnnotList[-1][1])
                iAnnotList[-1][2] = int(iAnnotList[-1][2])

    #Contig, Name, AbsStart, AbsEnd, CoordStart, CoordEnd, Identity, ..., relability
    finalDec = []
    with open(args.SDout, 'r') as f:
        csv_reader = csv.reader(f, delimiter='\t')
        for line in csv_reader:
            finalDec.append([line[0], line[1], cen_start + int(line[2]), cen_start + int(line[3]), int(line[2]), int(line[3]), line[4], line[-1]])

    #Contig, Name, Start, End, Name, Score, Starnd, Length, Gap, IsSamePos, Contig, Name, AbsStart, AbsEnd, Start, End, Identity, Relability
    merge_info = []
    curI = 0
    curJ = 0
    while curI < len(iAnnotList) or curJ < len(finalDec):
        if curI < len(iAnnotList) and (curJ == len(finalDec) or iAnnotList[curI][2] < finalDec[curJ][2]):
            merge_info.append(iAnnotList[curI] + ["FALSE"] + ["","","","","","","",""])
            curI += 1
        elif curJ < len(finalDec) and (curI == len(iAnnotList) or finalDec[curJ][3] < iAnnotList[curI][1]):
            merge_info.append(["","","","","","","","",""] + ["FALSE"] + finalDec[curJ])
            curJ += 1
        else:
            if iAnnotList[curI][1] == finalDec[curJ][2]:
                merge_info.append(iAnnotList[curI] + ["TRUE"] + finalDec[curJ])
            else:
                merge_info.append(iAnnotList[curI] + ["FALSE"] + finalDec[curJ])
            curJ += 1
            curI += 1

    with open(args.o, 'w') as fw:
        writer = csv.writer(fw)
        writer.writerow(["Contig", "Name", "Start", "End", "Name", "Score", "Starnd", "Length", "Gap", "IsSamePos", "Contig", "Name", "AbsStart", "AbsEnd", "Start", "End", "Identity", "Relability"])
        for row in merge_info:
            writer.writerow(row)


if __name__ == "__main__":
    main()