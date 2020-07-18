# (c) 2020 by Authors
# This file is a part of the SD program.
# see LICENSE file

from Bio import SeqIO


def read_bio_seq(filename):
    seqs = read_bio_seqs(filename)
    return str(list(seqs.values())[0])


def read_bio_seqs(filename):
    form = filename.split('.')[-1]
    if form == 'fa' or form == 'fna':
        form = 'fasta'
    elif form == 'fq':
        form = 'fastq'
    seqs = SeqIO.parse(filename, format=form)
    seqs = {seq.id: str(seq.seq) for seq in seqs}
    return seqs


trans = str.maketrans('ATGCatgc-', 'TACGtacg-')


def RC(s):
    return s.translate(trans)[::-1]


def write_bio_seqs(filename, seqs, width=60):
    with open(filename, 'w') as f:
        for seq_id, seq in seqs.items():
            print(f'>{seq_id}', file=f)
            if width is None:
                print(seq, file=f)
                continue
            start = 0
            seq_len = len(seq)
            while seq_len - start >= width:
                print(seq[start:start+width], file=f)
                start += width
            print(seq[start:], file=f)
