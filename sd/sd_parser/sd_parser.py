# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

import logging

import pandas as pd

from sd.monomers.monomer_db import MonomerDB
from sd.monomers.monostring_set import MonoStringSet
from sd.utils.bio import read_bio_seqs

logger = logging.getLogger("SD.sd_parser.sd_parser")


class SD_Report:
    def __init__(self, sd_report_fn, monomers_fn, sequences_fn,
                 mode=None, hpc=True, cluster=True):
        logger.info('Reading SD Report')

        if mode is not None:
            assert mode in ['ont', 'hifi', 'assembly']

        logger.info(f'Mode is {mode}')

        logger.info(f'    sd_report_fn = {sd_report_fn}')
        logger.info(f'    monomers_fn  = {monomers_fn}')
        logger.info(f'    sequences_fn = {sequences_fn}')
        monomer_db = MonomerDB.from_fasta_file(monomers_fn, cluster=cluster)

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
                                         mode=mode)

        logger.info('Creating monostrings dict')
        logger.info('Finished creating monostrings dict')

        self.monomer_db = monomer_db
        self.report = report
        self.monostring_set = monostring_set
        self.has_hpc = has_hpc
