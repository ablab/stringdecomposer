import argparse
from collections import namedtuple
import os
import sys

import edlib

SCRIPT_FN = os.path.realpath(__file__)
SCRIPT_DIR = os.path.dirname(SCRIPT_FN)
sys.path.insert(0, os.path.join(SCRIPT_DIR, os.pardir, os.pardir))

from sd.standard_logger import get_logger
from sd.sd_parser.sd_parser import run_SD, SD_Report
from sd.utils.bio import RC
from sd.utils.git import get_git_revision_short_hash
from sd.utils.os_utils import expandpath, smart_makedirs


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--target", help="Target sequences to encode",
                        required=True)
    parser.add_argument("-m", "--monomers", help="Monomers", required=True)
    parser.add_argument("--sd", help="SD report")
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("-t", "--threads", type=int, default=16)
    params = parser.parse_args()

    params.target = expandpath(params.target)
    params.monomers = expandpath(params.monomers)
    params.outdir = expandpath(params.outdir)
    if not os.path.isfile(params.target):
        logger.error(f'File does not exists --target == {params.target}')
        sys.exit(1)
    if not os.path.isfile(params.monomers):
        logger.error(f'File does not exists --monomers == {params.monomers}')
        sys.exit(1)
    return params


class EncodedChar:
    # char in seq can be '-', so it is not a nucleotide
    HEADER = ['mono_id', 'mono_is_reversed', 'monomer_pos', 'monomer_nucl',
              'seq_pos', 'seq_char', 'enc_seq_char']

    def __init__(self, mono_id, mono_is_reversed,
                 monomer_pos, monomer_nucl,
                 seq_pos, seq_char, enc_seq_char):
        self.mono_id = mono_id
        self.mono_is_reversed = mono_is_reversed
        self.monomer_pos = monomer_pos
        self.monomer_nucl = monomer_nucl
        self.seq_pos = seq_pos
        self.seq_char = seq_char
        self.enc_seq_char = enc_seq_char

    def output(self, opened_writable_file, sep='\t', print_header=True):
        if print_header:
            print(sep.join(self.HEADER), file=opened_writable_file)
        out_lst = [self.mono_id, self.mono_is_reversed, self.monomer_pos,
                   self.monomer_nucl, self.seq_pos, self.seq_char,
                   self.enc_seq_char]
        out_lst = [str(x) for x in out_lst]
        out_str = sep.join(out_lst)
        print(out_str, file=opened_writable_file)


class EncodedMonoString:
    def __init__(self, monostring, lst_encoded_chars):
        self.monostring = monostring
        self.lst_encoded_chars = lst_encoded_chars

    def output(self, outfile, sep='\t', print_header=True):
        with open(outfile, 'a') as f:
            if print_header:
                print(sep.join(EncodedChar.HEADER), file=f)
            for enc_char in self.lst_encoded_chars:
                enc_char.output(f, sep=sep, print_header=False)


class TandemEncoding:
    PositionTriple = namedtuple('PositionTriple', 'MonoIndex coord base')

    def __init__(self, monomer_db):
        self.monomer_db = monomer_db
        self.encoding = {}
        self.decoding = []
        i = 0
        bases = list("ACGT-")
        for monomer in monomer_db.monomers:
            mono_index = monomer.mono_index
            length = len(monomer.seq)
            for j in range(length):
                for nucl in bases:
                    triple = self.PositionTriple(mono_index, j, nucl)
                    self.encoding[triple] = i
                    self.decoding.append(triple)
                    i += 1

    def encode_char(self, mono_index, mono_is_reversed,
                    monomer_pos, monomer_nucl,
                    seq_pos, seq_char):
        enc_seq_nucl = self.encoding[self.PositionTriple(mono_index,
                                                         monomer_pos,
                                                         monomer_nucl)]
        mono_id = self.monomer_db.index2id[mono_index][0]
        return EncodedChar(mono_id, mono_is_reversed,
                           monomer_pos, monomer_nucl,
                           seq_pos, seq_char,
                           enc_seq_nucl)

    def decode_char(self, i):
        return self.decoding[i]

    def encode_monostring(self, monostring, trim=3):
        lst_encoded_chars = []
        for i, mi in enumerate(monostring.monoinstances):
            if not mi.is_reliable():
                continue
            monomer = mi.monomer
            mono_index = monomer.mono_index
            m_cons = self.monomer_db.get_seqs_by_index(mono_index)
            m_cons = m_cons[trim:-trim]
            if mi.is_reverse():  # mi.strand == Strand.REVERSE
                m_cons = RC(m_cons)
            nucl_segment = mi.nucl_segment

            alignment = edlib.align(query=m_cons, target=nucl_segment,
                                    mode='HW',
                                    task='path')

            cons_align, _, segm_align = \
                edlib.getNiceAlignment(alignment,
                                       m_cons,
                                       nucl_segment).values()
            if mi.is_reverse():
                monomer_pos = trim + len(m_cons) - 1
            else:
                monomer_pos = trim

            if len(alignment['locations']) > 1:
                sts = [loc[0] for loc in alignment['locations']]
                sts = list(set(sts))
                if len(sts) > 1:
                    print(alignment['locations'])
                assert len(sts) == 1  # TODO investigate if that is false
            assert len(alignment['locations']) > 0

            mi_st_alignment = alignment['locations'][0][0]
            seq_pos = mi.st + mi_st_alignment

            for cons_char, segm_char in zip(cons_align, segm_align):
                if cons_char == '-':
                    seq_pos += 1
                    continue

                assert (mi.is_reverse() and \
                        m_cons[-(monomer_pos - trim)-1] == cons_char) or \
                       (not mi.is_reverse() and \
                        m_cons[monomer_pos - trim] == cons_char)

                lst_encoded_chars.append(self.encode_char(
                    mono_index=mono_index,
                    mono_is_reversed=mi.is_reverse(),
                    monomer_pos=monomer_pos,
                    monomer_nucl=cons_char,
                    seq_pos=seq_pos,
                    seq_char=segm_char)
                )
                if mi.is_reverse():
                    monomer_pos -= 1
                else:
                    monomer_pos += 1

                if segm_char != '-':
                    assert monostring.nucl_sequence[seq_pos] == segm_char
                    seq_pos += 1

            assert (mi.is_reverse() and monomer_pos == trim - 1) or \
                   (not mi.is_reverse() and monomer_pos == trim + len(m_cons))

        encoded_monostring = \
            EncodedMonoString(monostring=monostring,
                              lst_encoded_chars=lst_encoded_chars)
        return encoded_monostring

    def output(self, outdir, sep='\t'):
        outfile = os.path.join(outdir, 'tandem_encoding.tsv')
        logger.info(f'Exporting tandem encoding into {outfile}')
        with open(outfile, 'w') as f:
            print(sep.join(['monomer', 'coord', 'base', 'enc']), file=f)
            for pos_triple, enc in self.encoding.items():
                monoid = self.monomer_db.index2id[pos_triple.MonoIndex][0]
                print(sep.join([monoid,
                                str(pos_triple.coord),
                                pos_triple.base,
                                str(enc)]),
                      file=f)


def tandem_encode(target_fn, monomers_fn, outdir, n_threads,
                  sd_report_fn=None):
    smart_makedirs(outdir)
    if sd_report_fn is None:
        sd_report_fn = run_SD(sequences_fn=target_fn,
                              monomers_fn=monomers_fn,
                              outdir=os.path.join(outdir, 'SD_report'),
                              n_threads=n_threads)

    sd_report = SD_Report(sd_report_fn=sd_report_fn,
                          monomers_fn=monomers_fn,
                          sequences_fn=target_fn,
                          attempt_reversing=False)

    monostring_set = sd_report.monostring_set
    monomer_db = monostring_set.monomer_db

    tandem_encoding = TandemEncoding(monomer_db)
    tandem_encoding.output(outdir)

    encoded_outfn = os.path.join(outdir, 'encoded_monostrings.tsv')
    logger.info(f'Exporting encoded monostrings in {encoded_outfn}')
    print_header = True
    for monostring in monostring_set.values():
        encoded_monostring = tandem_encoding.encode_monostring(monostring)
        encoded_monostring.output(encoded_outfn, print_header=print_header)
        print_header = False


def main():
    params = parse_args()
    smart_makedirs(params.outdir)
    logfn = os.path.join(params.outdir, 'tandem_encoder.log')
    global logger
    logger = get_logger(logfn,
                        logger_name='SD')
    logger.info(f'TandemEncoding with {SCRIPT_FN} started')
    logger.info('cmd: {}'.format(sys.argv))
    logger.info('git hash: {}'.format(get_git_revision_short_hash()))

    tandem_encode(target_fn=params.target,
                  monomers_fn=params.monomers,
                  outdir=params.outdir,
                  n_threads=params.threads,
                  sd_report_fn=params.sd)


if __name__ == "__main__":
    main()
