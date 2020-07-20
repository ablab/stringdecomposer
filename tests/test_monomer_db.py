import os
import unittest

from sd.monomers.monomer_db import MonomerDB


this_dirname = os.path.dirname(os.path.realpath(__file__))
monomers_fn = os.path.join(this_dirname, os.path.pardir,
                           'test_data', 'DXZ1_star_monomers.fa')


class TestMonomerDB(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.monomer_db = MonomerDB.from_fasta_file(fn=monomers_fn)
        super(TestMonomerDB, self).__init__(*args, **kwargs)

    def test_size(self):
        self.assertEqual(self.monomer_db.get_size(), 12)
