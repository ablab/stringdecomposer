# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

import logging
import os
from subprocess import check_call

import pandas as pd

from sd.monomers.monomer_db import MonomerDB
from sd.monomers.monostring_set import MonoStringSet
from sd.utils.bio import read_bio_seqs
from sd.utils.os_utils import smart_makedirs

logger = logging.getLogger("SD.sd_parser.sd_parser")


class SD_Report:
    def __init__(self, sd_report_fn, monomers_fn, sequences_fn,
                 mode=None, tocluster=False,
                 attempt_reversing=True):
        logger.info('Reading SD Report')

        if mode is not None:
            assert mode in ['ont', 'hifi', 'assembly']

        logger.info('Mode is {}'.format(mode))

        logger.info('    sd_report_fn = {}'.format(sd_report_fn))
        logger.info('    monomers_fn  = {}'.format(monomers_fn))
        logger.info('    sequences_fn = {}'.format(sequences_fn))
        monomer_db = MonomerDB.from_fasta_file(monomers_fn,
                                               tocluster=tocluster)

        logger.info('Reading SD Report from csv')
        report = pd.read_csv(sd_report_fn, sep='\t',
                             header=None)
        _, ncols = report.shape
        if ncols == 12:
            # report has HPC
            has_hpc = True
            names = ['s_id', 'monomer',
                     's_st', 's_en',
                     'identity',
                     'sec_monomer', 'sec_identity',
                     'homo_monomer', 'homo_identity',
                     'sec_homo_monomer', 'sec_homo_identity',
                     'reliability'
                     ]
        elif ncols == 8:
            has_hpc = False
            names = ['s_id', 'monomer',
                     's_st', 's_en',
                     'identity',
                     'sec_monomer', 'sec_identity',
                     'reliability'
                     ]
        else:
            assert False
        report.columns = names
        logger.info('Finished reading SD Report from csv')

        logger.info('Reading sequences')
        sequences = read_bio_seqs(sequences_fn)
        logger.info('Finished reading sequences')
        monostring_set = \
            MonoStringSet.from_sd_report(report=report,
                                         sequences=sequences,
                                         monomer_db=monomer_db,
                                         mode=mode,
                                         attempt_reversing=attempt_reversing)

        logger.info('Creating monostrings dict')
        logger.info('Finished creating monostrings dict')

        self.monomer_db = monomer_db
        self.report = report
        self.monostring_set = monostring_set
        self.has_hpc = has_hpc


def run_SD(sequences_fn, monomers_fn, outdir='.',
           outfn='final_decomposition.tsv',
           n_threads=4):
    logger.info(f'Running SD on')
    logger.info(f'\tSequences = {sequences_fn}')
    logger.info(f'\tMonomers = {monomers_fn}')
    logger.info(f'\tOutdir = {outdir}')
    smart_makedirs(outdir)
    outfn = os.path.join(outdir, outfn)
    if os.path.isfile(outfn):
        logger.info(f'File {outfn} exists. Reusing')
        return outfn

    script_dir = os.path.dirname(os.path.realpath(__file__))
    sd_bin = os.path.join(script_dir, os.pardir, 'run_decomposer.py')
    cmd = f'{sd_bin} {sequences_fn} {monomers_fn} ' + \
          f'-t {n_threads} -o {outfn}'
    logger.info(cmd)
    cmd = cmd.split(' ')
    check_call(cmd)
    return outfn
