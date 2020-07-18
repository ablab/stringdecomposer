# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

from collections import Counter
import logging

import numpy as np

from monomers.monostring import MonoString
from utils.bio import compress_homopolymer

logger = logging.getLogger("SD.monomers.monostring_set")


class MonoStringSet:
    def __init__(self, monostrings, monostrings_filt_out, monomer_db, mode):
        self.monostrings = monostrings
        self.monostrings_filt_out = monostrings_filt_out
        self.monomer_db = monomer_db
        self.mode = mode
        self.hybrids_corrected = False
        self.get_stats()

    @classmethod
    def from_sd_report(cls, report, sequences, monomer_db,
                       mode=None, filter_lowq=False):
        def get_raw_monostrings(report=report,
                                sequences=sequences,
                                monomer_db=monomer_db):
            monostrings = {}
            for seq_id, seq_record in report.groupby('s_id'):
                seq_record = seq_record.sort_values(by=['s_st'])
                nucl_sequence = sequences[seq_id]
                monostring = \
                    MonoString.from_sd_record(seq_id=seq_id,
                                              monomer_db=monomer_db,
                                              sd_record=seq_record,
                                              nucl_sequence=nucl_sequence)
                monostrings[seq_id] = monostring
            return monostrings

        def filter_monostrings(monostrings,
                               max_lowercase=None,
                               max_unreliable=None):
            if max_lowercase is None or max_unreliable is None:
                if mode == 'assembly':
                    max_lowercase = 1
                    max_unreliable = 1
                else:
                    max_lowercase = 0.1
                    max_unreliable = 0.4
            good_monostrings, bad_monostrings = {}, {}
            for seq_id, ms in monostrings.items():
                if ms.get_perc_lowercase() < max_lowercase and \
                   ms.get_perc_unreliable() < max_unreliable:
                    good_monostrings[seq_id] = ms
                else:
                    bad_monostrings[seq_id] = ms
            return good_monostrings, bad_monostrings

        if mode is not None:
            assert mode in ['ont', 'hifi', 'assembly']
        unfiltered_monostrings = get_raw_monostrings()
        if filter_lowq:
            assert mode is not None
            monostrings, monostrings_filt_out = \
                filter_monostrings(unfiltered_monostrings)
        else:
            monostrings = unfiltered_monostrings
            monostrings_filt_out = {}
        monostring_set = cls(monostrings=monostrings,
                             monostrings_filt_out=monostrings_filt_out,
                             monomer_db=monomer_db,
                             mode=mode)
        return monostring_set

    def get_nucl_seq(self, seq_id):
        if seq_id in self.monostrings:
            return self.monostrings[seq_id].nucl_sequence
        else:
            assert seq_id in self.monostrings_filt_out
            return self.monostrings_filt_out[seq_id].nucl_sequence

    def __getitem__(self, sub):
        return self.monostrings[sub]

    def __len__(self):
        return len(self.monostrings)

    def classify_monomerinstances(self, only_reliable=True):
        monoindexes = self.monomer_db.get_monoindexes()
        all_monomerinstances_dict = \
            {monoindex: [] for monoindex in monoindexes}
        for ms in self.monostrings.values():
            mi_dict = \
                ms.classify_monomerinstances(only_reliable=only_reliable)
            for monoindex in monoindexes:
                all_monomerinstances_dict[monoindex] += mi_dict[monoindex]
        return all_monomerinstances_dict

    def get_monomerinstances_by_monoindex(self, mono_index,
                                          only_reliable=True):
        all_monomerinstances_dict = \
            self.classify_monomerinstances(only_reliable=only_reliable)
        return all_monomerinstances_dict[mono_index]

    def get_stats(self, return_stats=False):
        def get_ngap_symbols(monostrings, compr_hmp=False):
            cnt = 0
            for monostring in monostrings.values():
                string = monostring.raw_monostring
                if compr_hmp:
                    string = compress_homopolymer(string, return_list=True)
                char_counter = Counter(string)
                cnt += char_counter[monostring.gap_symb]
            return cnt

        stats = {}
        monostrings = self.monostrings
        monostrings_filt_out = self.monostrings_filt_out
        monostrings_lens = [len(monostr) for monostr in monostrings.values()]
        stats['nmonostrings'] = len(monostrings_lens)
        stats['nfiltered_out'] = len(monostrings_filt_out)
        stats['min_len'] = np.min(monostrings_lens)
        stats['max_len'] = np.max(monostrings_lens)
        stats['mean_len'] = np.mean(monostrings_lens)
        stats['tot_len'] = np.sum(monostrings_lens)
        stats['ngaps'] = get_ngap_symbols(monostrings)
        stats['pgaps'] = stats['ngaps'] / stats['tot_len']
        stats['ngap_runs'] = get_ngap_symbols(monostrings, compr_hmp=True)

        logger.info(f'# monostrings: {stats["nmonostrings"]}')
        logger.info(f'# filtered monostrings: {stats["nfiltered_out"]}')
        logger.info(f'Min length = {stats["min_len"]}')
        logger.info(f'Mean length = {stats["mean_len"]}')
        logger.info(f'Max length = {stats["max_len"]}')
        logger.info(f'Total length = {stats["tot_len"]}')

        logger.info(f'#(%) Gap symbols = {stats["ngaps"]} ({stats["pgaps"]})')
        logger.info(f'#Gap runs = {stats["ngap_runs"]}')
        if return_stats:
            return stats

    def keys(self):
        return self.monostrings.keys()

    def values(self):
        return self.monostrings.values()

    def items(self):
        return self.monostrings.items()
