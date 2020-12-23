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


class EncodedMonoString:
    def __init__(self, monostring, lst_encoding):
        self.monostring = monostring
        self.lst_encoding = lst_encoding

    def output(self, outfile):
        with open(outfile, 'a') as f:
            print(f'>{self.monostring.seq_id}', file=f)
            print(' '.join(str(x) for x in self.lst_encoding), file=f)


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

    def encode_char(self, mono_index, pos, nucl):
        return self.encoding[self.PositionTriple(mono_index, pos, nucl)]

    def decode_char(self, i):
        return self.decoding[i]

    def encode_monostring(self, monostring, trim=3):
        lst_encoding = []
        for i, mi in enumerate(monostring.monoinstances):
            if not mi.is_reliable():
                continue
            monomer = mi.monomer
            mono_index = monomer.mono_index
            m_cons = self.monomer_db.get_seqs_by_index(mono_index)
            m_cons = m_cons[trim:-trim]
            nucl_segment = mi.nucl_segment
            alignment = edlib.align(m_cons, nucl_segment,
                                    mode='HW',
                                    task='path')
            qa, _, ta = edlib.getNiceAlignment(alignment,
                                               m_cons,
                                               nucl_segment).values()
            j = 0
            for qc, tc in zip(qa, ta):
                if qc == '-':
                    continue
                lst_encoding.append(self.encode_char(mono_index=mono_index,
                                                     pos=j,
                                                     nucl=qc))
                j += 1

        encoded_monostring = EncodedMonoString(monostring=monostring,
                                               lst_encoding=lst_encoding)
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
                          sequences_fn=target_fn)
    monostring_set = sd_report.monostring_set
    monomer_db = monostring_set.monomer_db

    tandem_encoding = TandemEncoding(monomer_db)
    tandem_encoding.output(outdir)

    encoded_outfn = os.path.join(outdir, 'encoded_monostrings.extfa')
    logger.info(f'Exporting encoded monostrings in {encoded_outfn}')
    for monostring in monostring_set.values():
        encoded_monostring = tandem_encoding.encode_monostring(monostring)
        encoded_monostring.output(encoded_outfn)


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
