# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

from collections import defaultdict
import logging

from sd.cluster_sequences import cluster_sequences
from sd.utils.bio import read_bio_seqs
from sd.utils.os_utils import expandpath

logger = logging.getLogger("SD.monomers.monomer_db")


class Monomer:
    def __init__(self, monomer_id, mono_index, seq):
        self.monomer_id = monomer_id
        self.mono_index = mono_index
        self.seq = seq

    def __repr__(self):
        return f'monomer_id={self.monomer_id}, '\
               f'mono_index={self.mono_index}, '\
               f'seq={self.seq}'


class MonomerDB:
    def __init__(self, id2index, index2id, monomers, id2list_coord, clustered):
        self.id2index = id2index
        self.index2id = index2id
        self.monomers = monomers
        self.id2list_coord = id2list_coord
        self.clustered = clustered

    def __repr__(self):
        return f'size={self.get_size()}, ids={self.get_ids()}'

    @classmethod
    def from_fasta_file(cls, fn, cluster_max_ident=0.95, tocluster=False):
        fn = expandpath(fn)
        logger.info('Creating Monomer DataBase from {}'.format(fn))
        raw_monomers = read_bio_seqs(fn)
        logger.info('Clustering monomers.'
                    'Identity thresh {}'.format(cluster_max_ident))
        if tocluster:
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
                logger.debug('Monomer: index = {} id = {}'.format(i,
                                                                  monomer_id))
                logger.debug('         monomer seq = {}'.format(monomer_seq))

        monomer_db = cls(id2index=id2index,
                         index2id=index2id,
                         monomers=monomers,
                         id2list_coord=id2list_coord,
                         clustered=tocluster)

        logger.info('Finished Creating Monomer DataBase')
        return monomer_db

    def get_ids(self):
        return self.id2index.keys()

    def get_monoindexes(self):
        return self.index2id.keys()

    def get_monomer_by_id(self, mono_id):
        return self.monomers[self.id2index[mono_id]]

    def get_monomers_by_index(self, mono_index):
        monomers = []
        for mono_id in self.index2id[mono_index]:
            list_coord = self.id2list_coord[mono_id]
            monomers.append(self.monomers[list_coord])
        if not self.clustered:
            assert len(monomers) == 1
            monomer = monomers[0]
            return monomer
        return monomers

    def get_monomers_dict(self):
        return {monomer.monomer_id: monomer.seq
                for monomer in self.monomers}

    def get_seq_by_id(self, monomer_id):
        mono_index = self.id2index[monomer_id]
        return self.monomers[mono_index].seq

    def get_seqs_by_index(self, mono_index):
        assert 0 <= mono_index <= max(self.index2id)
        if not self.clustered:
            monomer = self.get_monomers_by_index(mono_index)
            return monomer.seq
        seqs = []
        for monomer in self.get_monomers_by_index(mono_index):
            seqs.append(monomer.seq)
        return seqs

    def get_size(self):
        return 1 + max(self.index2id)
