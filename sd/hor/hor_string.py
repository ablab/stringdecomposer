# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

import itertools
import logging

from sd.hor.hor_monosequence import HORMonoSequence

logger = logging.getLogger("SD.hor.hor_string")


class HORInstance:
    def __init__(self, hor_monosequence,
                 monoinstances,
                 identity,
                 st, en,
                 nucl_segment,
                 seq_len,
                 is_NM):
        assert en - st == len(nucl_segment)
        self.hor_monosequence = hor_monosequence
        self.monoinstances = monoinstances
        self.identity = identity
        self.st = st
        self.en = en
        self.nucl_segment = nucl_segment
        self.seq_len = seq_len
        self.is_NM = is_NM

        for minst in monoinstances:
            assert minst.seq_len == seq_len
            assert self.nucl_segment[minst.st-self.st:minst.en-self.st] == \
                minst.nucl_segment

    def get_sequence_mono_ids(self):
        if self.is_NM:
            return tuple('?')
        return self.hor_monosequence.mono_ids

    def get_mono_indexes(self):
        if self.is_NM:
            return None
        return self.hor_monosequence.get_mono_indexes()


class HORString:
    NM = 'NM'

    def __init__(self, seq_id, hor_instances, monostring):
        self.seq_id = seq_id
        self.hor_instances = hor_instances
        self.monostring = monostring

    @classmethod
    def from_hor_record(cls, seq_id, record, monostring):
        record.s_en += 1
        record.hor_monosequence = \
            record.hor_monosequence.map(lambda x: x.split(','))

        record_nucl_len = len(monostring.nucl_sequence)
        if monostring.is_reversed:
            record.s_st, record.s_en = \
                record_nucl_len - record.s_en, record_nucl_len - record.s_st
            record = record[::-1]
            record.hor_monosequence = \
                record.hor_monosequence.map(lambda x: x[::-1])

        starts = list(record.s_st)
        ends = list(record.s_en)
        identities = [ident / 100
                      for ident in record.identity]
        hor_mss = record.hor_monosequence

        hor_instances = []
        i = 0
        for st, en, ident, hor_ms in zip(starts, ends, identities, hor_mss):
            while i < len(monostring) and st != monostring.monoinstances[i].st:
                i += 1
            j = i
            while j < len(monostring) and en != monostring.monoinstances[j].en:
                j += 1
            j += 1
            monoinstances = monostring.monoinstances[i:j]
            is_NM = hor_ms[0] == cls.NM
            if is_NM:
                hor_ms = None
            else:
                mis_ids = tuple(mi.get_monoid() for mi in monoinstances)
                # for short_mono_id, mi in zip(hor_ms, mis_ids):
                #     assert short_mono_id == mi[0]
                hor_ms = HORMonoSequence(mono_ids=mis_ids,
                                         monomer_db=monostring.monomer_db)
            nucl_segment = monostring.get_nucl_segment(st, en)
            hor_instance = HORInstance(hor_monosequence=hor_ms,
                                       monoinstances=monoinstances,
                                       identity=ident,
                                       st=st, en=en,
                                       nucl_segment=nucl_segment,
                                       seq_len=len(monostring.nucl_sequence),
                                       is_NM=is_NM)
            hor_instances.append(hor_instance)
            i = j

        hor_string = cls(seq_id=seq_id,
                         hor_instances=hor_instances,
                         monostring=monostring)
        return hor_string

    def __len__(self):
        return self.hor_instances.__len__()

    def __getitem__(self, sub):
        if isinstance(sub, slice):
            sublist = self.hor_instances[sub.start:sub.stop:sub.step]
            sublist = [hor_instance.get_mono_indexes()
                       for hor_instance in sublist]
            return sublist
        return self.hor_instances[sub].get_mono_indexes()

    def get_compact_form(self):
        raw_hor_string = []
        for hor_instance in self.hor_instances:
            mono_ids = hor_instance.get_sequence_mono_ids()
            if len(mono_ids) > 1:
                raw_hor_string.append(f'{mono_ids[0][0]}-{mono_ids[-1][0]}')
            else:
                raw_hor_string.append(f'{mono_ids[0][0]}')
        raw_hor_string = itertools.groupby(raw_hor_string)
        raw_hor_string = [(e, len(list(g))) for e, g in raw_hor_string]
        return raw_hor_string

    def get_nucl_segment(self, i):
        return self.hor_instances[i].nucl_segment

    def is_reversed(self):
        return self.monostring.is_reversed

    def get_hor_instance_coords(self, i, original=False):
        assert 0 <= i < len(self.hor_instances)
        st = self.hor_instances[i].st
        en = self.hor_instances[i].en
        if original and self.is_reversed():
            seq_len = self.monostring.get_seq_len()
            st, en = seq_len - st, seq_len - en
        return st, en
