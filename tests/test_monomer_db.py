import os
import unittest

from sd.monomers.monomer_db import MonomerDB
from sd.utils.various import fst_iterable


this_dirname = os.path.dirname(os.path.realpath(__file__))
monomers_fn = os.path.join(this_dirname, os.path.pardir,
                           'test_data', 'DXZ1_star_monomers.fa')


class TestMonomerDB(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.monomer_db = MonomerDB.from_fasta_file(fn=monomers_fn)
        self.mono_id = "A_0_DXZ1*_doubled/1978_2147/R"
        self.mono_index = 0
        self.mono_seq = "TCCGTTTAGCTTTTAGGTGAAGATTATCCCGTTTCCAACGAAACCTTCAAAG"\
                        "AGGTCCAAATATCCCCTTGCGGATCCCACAGAAAGAGTGTTTCGAAACTGCT"\
                        "GTTTCAAAGGAATCTTCAACTCTGTGAGTTGAATGCAATCATCACAAAGAAG"\
                        "TTTCTGACAATGCT"
        self.db_size = 12
        super(TestMonomerDB, self).__init__(*args, **kwargs)

    def test_get_ids(self):
        ids = self.monomer_db.get_ids()
        self.assertEqual(len(ids), self.db_size)
        self.assertEqual(fst_iterable(ids), self.mono_id)

    def test_get_monoindexes(self):
        monoindexes = self.monomer_db.get_monoindexes()
        self.assertEqual(len(monoindexes),
                         self.db_size)
        self.assertEqual(list(monoindexes), list(range(self.db_size)))

    def test_get_monomer_by_id(self):
        monomer = self.monomer_db.get_monomer_by_id(self.mono_id)
        self.assertEqual(monomer.monomer_id, self.mono_id)

    def test_get_monomers_dict(self):
        self.monomer_db.get_monomers_dict()

    def test_get_seq_by_id(self):
        seq = self.monomer_db.get_seq_by_id(self.mono_id)
        self.assertEqual(seq, self.mono_seq)

    def test_get_seq_by_index(self):
        seq = self.monomer_db.get_seqs_by_index(self.mono_index)
        self.assertEqual(seq, self.mono_seq)

    def test_size(self):
        self.assertEqual(self.monomer_db.get_size(), self.db_size)
