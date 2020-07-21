# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

import logging

from sd.monomers.monomer_db import MonomerDB

logger = logging.getLogger("SD.hor.hor")


# TODO: currently disabled class
class HOR:
    def __init__(self, mono_indexes, monomer_db):
        self.mono_indexes = mono_indexes
        self.monomer_db = monomer_db

    @classmethod
    def from_canonical_monomerdb(cls, canonical_fn,
                                 monomers_fn=None, monomer_db=None):
        assert (monomers_fn is None) != (monomer_db is None)
        if monomer_db is None:
            monomer_db = MonomerDB.from_fasta_file(fn=monomers_fn)

        with open(canonical_fn) as f:
            canonical_hor = f.readline().strip()
        canonical_hor = canonical_hor.split(',')
        logger.info('Canonical HOR: {}'.format(canonical_hor))
        mono_indexes = set(monomer_db.id2index[mono_id]
                           for mono_id in canonical_hor)
        hor = cls(mono_indexes=mono_indexes,
                  monomer_db=monomer_db)
        return hor

    def get_ids(self):
        return set(tuple(self.monomer_db.index2id[index])
                   for index in self.mono_indexes)

    def get_mono_indexes(self):
        return self.mono_indexes
