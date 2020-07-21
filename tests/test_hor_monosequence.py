import os
import unittest

from sd.utils.bio import read_bio_seq, RC
from sd.utils.various import fst_iterable

from sd.hor.hor import HOR
from sd.hor.hor_monosequence import HORMonoSequence
from sd.monomers.monomer_db import MonomerDB


this_dirname = os.path.dirname(os.path.realpath(__file__))
canonical_fn = os.path.join(this_dirname, os.path.pardir,
                            'test_data', 'canonical_X.txt')
monomers_fn = os.path.join(this_dirname, os.path.pardir,
                           'test_data', 'DXZ1_star_monomers.fa')
dxz1_star_fn = os.path.join(this_dirname, os.path.pardir,
                            'test_data', 'DXZ1_rc_star.fasta')


class TestHORMonoSequence(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.monomer_db = MonomerDB.from_fasta_file(fn=monomers_fn)
        dxz1_star = read_bio_seq(dxz1_star_fn)
        self.dxz1_star = dxz1_star[1977:] + dxz1_star[:1977]
        with open(canonical_fn) as f:
            canonical_hor = f.readline().strip()
        canonical_hor = canonical_hor.split(',')
        self.horms = HORMonoSequence(monomer_db=self.monomer_db,
                                     mono_ids=canonical_hor)
        super(TestHORMonoSequence, self).__init__(*args, **kwargs)

    def test_consensus(self):
        consensus = self.horms.get_consensus()
        self.assertEqual(self.dxz1_star, consensus)
