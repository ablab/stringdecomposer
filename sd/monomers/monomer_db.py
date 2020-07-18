# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

from collections import defaultdict
import logging

from utils.bio import read_bio_seqs
from utils.cluster_sequences import cluster_sequences
from utils.os_utils import expandpath

logger = logging.getLogger("SD.monomers.monomer_db")


class Monomer:
    def __init__(self, monomer_id, mono_index, seq):
        self.monomer_id = monomer_id
        self.mono_index = mono_index
        self.seq = seq

    def __repr__(self):
        return f'monomer_id={self.monomer_id}, mono_index={self.mono_index}, seq={self.seq}'


class MonomerDB:
    def __init__(self, id2index, index2id, monomers, id2list_coord):
        self.id2index = id2index
        self.index2id = index2id
        self.monomers = monomers
        self.id2list_coord = id2list_coord

    @classmethod
    def from_fasta_file(cls, fn, cluster_max_ident=0.95, cluster=True):
        fn = expandpath(fn)
        logger.info(f'Creating Monomer DataBase from {fn}')
        raw_monomers = read_bio_seqs(fn)
        logger.info(f'Clustering monomers.'
                    f'Identity thresh {cluster_max_ident}')
        if cluster:
            monomer_clusters, _ = \
                cluster_sequences(sequences=raw_monomers,
                                  max_ident=cluster_max_ident)
        else:
            # if no clustering, each sequence will be in its own cluster
            monomer_clusters = [[monomer_id] for monomer_id in raw_monomers]

        id2index = {}
        index2id = defaultdict(list)
        monomers = []
        id2list_coord = {}
        for i, cluster in enumerate(monomer_clusters):
            for monomer_id in cluster:
                monomer_seq = raw_monomers[monomer_id]
                monomer = Monomer(monomer_id=monomer_id,
                                  mono_index=i,
                                  seq=monomer_seq)
                monomers.append(monomer)
                id2index[monomer_id] = i
                index2id[i].append(monomer_id)
                id2list_coord[monomer_id] = len(monomers) - 1
                logger.debug(f'Monomer: index = {i} id = {monomer_id}')
                logger.debug(f'         monomer sequence = {monomer_seq}')

        monomer_db = cls(id2index=id2index,
                         index2id=index2id,
                         monomers=monomers,
                         id2list_coord=id2list_coord)

        logger.info(f'Finished Creating Monomer DataBase')
        return monomer_db

    def get_monomer_by_id(self, mono_id):
        return self.monomers[self.id2index[mono_id]]

    def get_monomers_by_index(self, mono_index):
        for mono_id in self.index2id[mono_index]:
            list_coord = self.id2list_coord[mono_id]
            yield self.monomers[list_coord]

    def get_seqs_by_index(self, mono_index):
        assert 0 <= mono_index <= max(self.index2id)
        for monomer in self.get_monomers_by_index(mono_index):
            yield monomer.seq

    def get_seq_by_id(self, monomer_id):
        mono_index = self.id2index[monomer_id]
        return self.monomers[mono_index].seq

    def get_monoindexes(self):
        return self.index2id.keys()

    def get_ids(self):
        return self.id2index.keys()

    def get_size(self):
        return 1 + max(self.index2id)

    def get_monomers_dict(self):
        return {monomer.monomer_id: monomer.seq
                for monomer in self.monomers}

    def extend_db(self, extra_monomers):
        index = max(self.index2id)
        for monomer_id, seq in extra_monomers.items():
            if monomer_id in self.id2index:
                continue
            index += 1
            monomer = Monomer(monomer_id=monomer_id,
                              mono_index=index,
                              seq=seq)
            self.id2index[monomer_id] = index
            self.index2id[index].append(monomer_id)
            self.monomers.append(monomer)
            self.id2list_coord[monomer_id] = len(self.monomers) - 1
