# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

import logging
import os
from subprocess import check_call

import pandas as pd

from sd.hor.hor_string_set import HORStringSet
from sd.sd_parser.sd_parser import SD_Report

from sd.utils.os_utils import expandpath, smart_makedirs

logger = logging.getLogger("SD.hor.hor_extraction_parser")


class HORExtractionReport:
    def __init__(self, hor_report_fn, sd_report_fn,
                 monomers_fn, sequences_fn,
                 tocluster=False):
        logger.info("Reading HOR extraction report")

        logger.info('    sd_report_fn  = {}'.format(sd_report_fn))
        logger.info('    hor_report_fn = {}'.format(hor_report_fn))
        logger.info('    monomers_fn   = {}'.format(monomers_fn))
        logger.info('    sequences_fn  = {}'.format(sequences_fn))

        sd_report = SD_Report(sd_report_fn=sd_report_fn,
                              monomers_fn=monomers_fn,
                              sequences_fn=sequences_fn)

        logger.info('Reading HOR extraction Report from csv')
        hor_report = pd.read_csv(hor_report_fn, sep='\t',
                                 header=None)
        names = ['s_id', 'hor_monosequence', 'n_monomers',
                 'identity', 's_st', 's_en', 'seq_len', 'gap_len']
        hor_report.columns = names

        logger.info('Finished reading HOR Extraction Report from csv')

        monomer_db = sd_report.monomer_db
        monostring_set = sd_report.monostring_set

        horstring_set = \
            HORStringSet.from_HOR_report(report=hor_report,
                                         monostring_set=monostring_set)
        self.monomer_db = monomer_db
        self.horstring_set = horstring_set
        self.hor_report = hor_report
        self.sd_report = sd_report


def run_hor_extractor(sequences_fn, sd_report_fn, monomers_fn, canonical_fn,
                      outdir,
                      outfn='hor_decomposition.tsv'):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    hor_extractor = os.path.join(script_dir, os.pardir,
                                  'scripts', 'extract_hors.py')
    hor_extractor = expandpath(hor_extractor)
    logger.info(f'Running HOR extraction with {hor_extractor} on')
    logger.info(f'\tSequences = {sequences_fn}')
    logger.info(f'\tMonomers = {monomers_fn}')
    logger.info(f'\tSD Report = {sd_report_fn}')
    logger.info(f'\tCanonical HORs = {canonical_fn}')
    logger.info(f'\tOutdir = {outdir}')

    smart_makedirs(outdir)
    outfn = os.path.join(outdir, outfn)
    if os.path.isfile(outfn):
        logger.info(f'File {outfn} exists. Reusing')
        return outfn

    cmd = f'python {hor_extractor} --canonical {canonical_fn} ' +\
          f'{sequences_fn} {monomers_fn} {sd_report_fn} {outfn}'
    logger.info(cmd)
    cmd = cmd.split(' ')
    check_call(cmd)
    return outfn
