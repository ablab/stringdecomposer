# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

import logging

logger = logging.getLogger("SD.hor.hor_monosequence")


class HORMonoSequence:
    def __init__(self, mono_ids, monomer_db):
        # TODO save an HOR
        # mono_ids is a list
        self.mono_ids = mono_ids
        self.monomer_db = monomer_db

    def get_mono_indexes(self):
        return tuple(self.monomer_db.id2index[mono_id]
                     for mono_id in self.mono_ids)

    def get_consensus(self):
        consensus = []
        for mono_id in self.mono_ids:
            consensus += self.monomer_db.get_seq_by_id(monomer_id=mono_id)
        consensus = ''.join(consensus)
        return consensus
