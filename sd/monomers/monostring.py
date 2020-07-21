# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

from enum import Enum
from itertools import count
import logging

import numpy as np

from sd.utils.bio import RC

logger = logging.getLogger("SD.monomers.monostring")


class Strand(Enum):
    FORWARD = '+'
    REVERSE = '-'

    @staticmethod
    def switch(strand):
        if strand == Strand.FORWARD:
            return Strand.REVERSE
        else:
            assert strand == Strand.REVERSE
            return Strand.FORWARD


class Reliability(Enum):
    RELIABLE = '+'
    UNRELIABLE = '?'


def assert_monostring_validity(monostring):
    string = monostring.raw_monostring
    monomer_db = monostring.monomer_db
    monomer_db_size = monomer_db.get_size()
    monoinsts = monostring.monoinstances
    for i, monoinstance in enumerate(monoinsts):
        mono_index = monoinstance.get_monoindex()
        if monoinstance.strand == Strand.REVERSE:
            mono_index += monomer_db_size
        if monoinstance.reliability == Reliability.RELIABLE:
            assert mono_index == string[i]

    nucl_sequence = monostring.nucl_sequence
    for mi in monoinsts:
        assert nucl_sequence[mi.st:mi.en] == mi.nucl_segment


class MonoInstance:
    def __init__(self, monomer, sec_monomer, strand, sec_strand,
                 seq_id, nucl_segment,
                 st, en, seq_len,
                 reliability,
                 identity, sec_identity):
        assert en - st == len(nucl_segment)
        self.monomer = monomer
        self.sec_monomer = sec_monomer
        self.strand = strand
        self.sec_strand = sec_strand
        self.seq_id = seq_id
        self.nucl_segment = nucl_segment
        self.st = st
        self.en = en
        self.seq_len = seq_len
        self.reliability = reliability
        self.identity = identity
        self.sec_identity = sec_identity
        self.is_reversed = False

    def get_monoid(self):
        return self.monomer.monomer_id

    def get_secmonoid(self):
        return self.sec_monomer.monomer_id

    def get_monoindex(self):
        return self.monomer.mono_index

    def get_secmonoindex(self):
        return self.sec_monomer.mono_index

    def get_ref_seq(self):
        return self.monomer.seq

    def is_reverse(self):
        return self.strand == Strand.REVERSE

    def is_forward(self):
        return self.strand == Strand.FORWARD

    def is_reliable(self):
        return self.reliability is Reliability.RELIABLE

    def reverse(self):
        self.nucl_segment = RC(self.nucl_segment)
        self.strand = Strand.switch(self.strand)
        self.sec_strand = Strand.switch(self.sec_strand)
        # [st; en)
        self.st, self.en = self.seq_len - self.en, self.seq_len - self.st
        self.is_reversed = not self.is_reversed


class MonoString:
    # monostring is stored as a tuple because
    # |monomer_db| can exceed |ascii|
    gap_symb = '?'
    reverse_symb = "'"
    none_monomer = 'None'

    def __init__(self, seq_id, monoinstances, raw_monostring, nucl_sequence,
                 monomer_db, is_reversed):
        self.seq_id = seq_id
        self.monoinstances = monoinstances
        self.raw_monostring = raw_monostring
        self.nucl_sequence = nucl_sequence
        self.monomer_db = monomer_db
        self.is_reversed = is_reversed
        self.corrections = {}  # dict pos (int) -> mono_index (int)
        assert_monostring_validity(self)

    @classmethod
    def from_sd_record(cls, seq_id, monomer_db, sd_record, nucl_sequence):
        def get_monoinstances(sd_record):
            def id2index_strand(monomer_id, monomer_db=monomer_db):
                if monomer_id == cls.none_monomer:
                    index = None
                    strand = Strand.FORWARD
                    return index, strand

                if monomer_id[-1] == cls.reverse_symb:
                    monomer_id = monomer_id[:-1]
                    strand = Strand.REVERSE
                else:
                    strand = Strand.FORWARD
                index = monomer_db.id2index[monomer_id]
                return index, strand

            def get_reliablities(sd_record, identities, sec_identities):
                reliabilities = []
                for raw_rel, ident, sec_ident in zip(sd_record.reliability,
                                                     identities,
                                                     sec_identities):
                    reliability = Reliability(raw_rel)
                    reliabilities.append(reliability)
                return reliabilities

            starts = list(sd_record.s_st)
            ends = [en + 1 for en in sd_record.s_en]

            ids = list(sd_record.monomer)
            indexes_strands = map(id2index_strand, ids)
            indexes, strands = zip(*indexes_strands)

            sec_ids = list(sd_record.sec_monomer)
            sec_indexes_strands = map(id2index_strand, sec_ids)
            sec_indexes, sec_strands = zip(*sec_indexes_strands)

            identities = [ident / 100
                          for ident in sd_record.identity]
            sec_identities = [ident / 100
                              for ident in sd_record.sec_identity]

            reliabilities = get_reliablities(sd_record=sd_record,
                                             identities=identities,
                                             sec_identities=sec_identities)

            monoinstances = []
            for i, st, en, \
                    rel, strand, sec_strand, \
                    mono_index, sec_mono_index, \
                    identity, sec_identity in \
                    zip(count(), starts, ends,
                        reliabilities, strands, sec_strands,
                        indexes, sec_indexes,
                        identities, sec_identities):
                monomer = monomer_db.monomers[mono_index]
                sec_monomer = None
                if sec_mono_index is not None:
                    sec_monomer = monomer_db.monomers[sec_mono_index]
                nucl_segment = nucl_sequence[st:en]

                monoinstance = MonoInstance(monomer=monomer,
                                            sec_monomer=sec_monomer,
                                            strand=strand,
                                            sec_strand=sec_strand,
                                            seq_id=seq_id,
                                            nucl_segment=nucl_segment,
                                            st=st,
                                            en=en,
                                            seq_len=len(nucl_sequence),
                                            reliability=rel,
                                            identity=identity,
                                            sec_identity=sec_identity)
                monoinstances.append(monoinstance)
            return monoinstances

        def reverse_if_needed(monoinstances, nucl_sequence,
                              max_reverse=0.5):
            is_reverse = [monoinstance.is_reverse()
                           for monoinstance in monoinstances
                           if monoinstance.is_reliable()]
            perc_reverse = np.mean(is_reverse)
            to_reverse = perc_reverse > max_reverse
            if to_reverse:
                # reverse monoinstance
                monoinstances.reverse()
                for monoinstance in monoinstances:
                    monoinstance.reverse()
                # reverse nucl_sequence
                nucl_sequence = RC(nucl_sequence)
            return monoinstances, nucl_sequence, to_reverse

        def get_string(monoinstance):
            string = []
            for monoinstance in monoinstances:
                mono_index = monoinstance.get_monoindex()
                if monoinstance.reliability == Reliability.RELIABLE:
                    if monoinstance.strand == Strand.FORWARD:
                        string.append(mono_index)
                    else:
                        assert monoinstance.strand == Strand.REVERSE
                        string.append(mono_index + monomer_db.get_size())
                else:
                    assert monoinstance.reliability == Reliability.UNRELIABLE
                    string.append(cls.gap_symb)
            string = tuple(string)
            return string

        monoinstances = get_monoinstances(sd_record=sd_record)

        monoinstances, nucl_sequence, is_reversed = \
            reverse_if_needed(monoinstances, nucl_sequence)
        string = get_string(monoinstances)

        monostring = cls(seq_id=seq_id,
                         monoinstances=monoinstances,
                         raw_monostring=string,
                         nucl_sequence=nucl_sequence,
                         monomer_db=monomer_db,
                         is_reversed=is_reversed)
        return monostring

    def __len__(self):
        return self.raw_monostring.__len__()

    def __getitem__(self, sub):
        if isinstance(sub, slice):
            sublist = self.raw_monostring[sub.start:sub.stop:sub.step]
            return sublist
        return self.raw_monostring[sub]

    def __setitem__(self, sub, item):
        # sub - position
        # item - mono_index
        assert self.raw_monostring[sub] != item
        self.corrections[sub] = (self.raw_monostring[sub], item)
        self.raw_monostring[sub] = item

    def classify_monomerinstances_by_monoindex(self, only_reliable=True):
        monoindexes = self.monomer_db.get_monoindexes()
        monomerinstances_dict = {monoindex: [] for monoindex in monoindexes}
        for mi in self.monoinstances:
            if (not only_reliable) or (only_reliable and mi.is_reliable()):
                monoindex = mi.get_monoindex()
                monomerinstances_dict[monoindex].append(mi)
        return monomerinstances_dict

    def get_monomerinstances_by_monoindex(self, mono_index,
                                          only_reliable=True):
        monomerinstances_dict = \
            self.classify_monomerinstances_by_monoindex(
                only_reliable=only_reliable)
        return monomerinstances_dict[mono_index]

    def get_nucl_segment(self, st, en):
        assert 0 <= st < en < len(self.nucl_sequence)
        return self.nucl_sequence[st:en]

    def get_identities(self):
        identities = []
        for mi in self.monoinstances:
            identities.append(mi.identity)
        return identities

    def get_perc_reliable(self):
        is_reliable = [monoinstance.is_reliable()
                       for monoinstance in self.monoinstances]
        perc_reliable = np.mean(is_reliable)
        return perc_reliable

    def get_perc_unreliable(self):
        return 1 - self.get_perc_reliable()

    def get_perc_forward_strand(self):
        is_forward = [monoinstance.is_forward()
                      for monoinstance in self.monoinstances
                      if monoinstance.is_reliable()]
        perc_forward = np.mean(is_forward)
        return perc_forward

    def get_perc_reverse_strand(self):
        return 1 - self.get_perc_forward_strand()

    def is_corrected(self):
        return len(self.corrections)
