# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

import logging

from sd.hor.hor_string import HORString

logger = logging.getLogger("SD.hor.hor_string_set")


class HORStringSet:
    def __init__(self, horstrings):
        self.horstrings = horstrings

    @classmethod
    def from_HOR_report(cls, report, monostring_set):
        horstrings = {}
        for seq_id, record in report.groupby('s_id'):
            record = record.sort_values(by=['s_st'])
            monostring = monostring_set[seq_id]
            horstring = \
                HORString.from_hor_record(seq_id=seq_id,
                                          record=record,
                                          monostring=monostring)
            horstrings[seq_id] = horstring

        horstring_set = HORStringSet(horstrings)
        return horstring_set

    def __getitem__(self, sub):
        return self.horstrings[sub]

    def __len__(self):
        return len(self.horstrings)

