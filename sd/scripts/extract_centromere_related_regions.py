import sys
import subprocess

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import edlib
import argparse

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

def save_fasta(filename, orfs):
    with open(filename, "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")

def edist_hw(lst):
    if len(str(lst[0])) == 0:
        return 100500
    if len(str(lst[1])) == 0:
        return 100500
    result = edlib.align(str(lst[0]), str(lst[1]), mode="HW", k = lst[2])
    if result["editDistance"] == -1:
        return 100500
    return result["editDistance"]


def find_as_borders_mono(seq, monomers, threshold):
    step = 10000
    start, end = len(seq), None
    for i in range(0, len(seq), step):
        min_ed = 100500
        for m in monomers:
            ed = edist_hw([m.seq, seq[i:i+step], threshold])
            min_ed = min(min_ed, ed)
            ed = edist_hw([m.seq.reverse_complement(), seq[i:i+step], threshold])
            min_ed = min(min_ed, ed)
        if min_ed < threshold:
            if start > i:
                start = i
            end = i + step
    return start, end


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Cuts reads, leaving only centromere-related part')
    parser.add_argument('-s', '--sequences', help='fasta-file with long reads sequences', required=True)
    parser.add_argument('-m', '--monomers', help='fasta-file with monomers', required=True)
    parser.add_argument('-o', '--out', help='path to output fasta-file', required=True)
    parser.add_argument('-d', '--edit-distance', help='threshold on edit distance for monomer alignment (default is 30)', default=30, required=False)
    args = parser.parse_args()

    reads_path = args.sequences
    mono_file = args.monomers
    monomers = load_fasta(mono_file)
    centromeric_reads = []

    total_len, centromeric_len, total_reads, centromeric_reads_num = 0, 0, 0, 0

    for read in SeqIO.parse(reads_path, "fasta"):
        start, end =  find_as_borders_mono(read.seq, monomers, args.edit_distance)
        total_len += len(read.seq)
        total_reads += 1
        if start == len(read.seq):
            print(read.id, " Filtered")
        else:
            centromeric_reads_num += 1
            centromeric_len += end - start + 1
            print(read.id, " Full read len=", len(read.seq), "Centromere start=", start, " Centromere end=", end, "Centromere part len=", end - start)
            centromeric_reads.append(make_record(read.seq[start:end + 1], read.id + ":" + str(start) + "-" + str(end), read.id + ":" + str(start) + "-" + str(end)))
    print("Total number of reads ", total_reads, " Total reads length ", total_len)
    print("Number of centromeric reads ", total_reads, " Centromeric reads length ", total_len)
    save_fasta(args.out, centromeric_reads)
    print("New seqeunces saved to ", args.out)