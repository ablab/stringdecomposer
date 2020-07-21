import os
import unittest

from sd.utils.various import fst_iterable

from sd.hor.hor import HOR


this_dirname = os.path.dirname(os.path.realpath(__file__))
canonical_fn = os.path.join(this_dirname, os.path.pardir,
                            'test_data', 'canonical_X.txt')
monomers_fn = os.path.join(this_dirname, os.path.pardir,
                           'test_data', 'DXZ1_star_monomers.fa')


class TestHOR(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.hor = HOR.from_canonical_monomerdb(canonical_fn=canonical_fn,
                                                monomers_fn=monomers_fn)
        super(TestHOR, self).__init__(*args, **kwargs)

    def test_get_ids(self):
        ids = self.hor.get_ids()
        self.assertEqual(fst_iterable(sorted(ids)),
                         ('A_0_DXZ1*_doubled/1978_2147/R',))

    def test_mono_indexes(self):
        self.assertEqual(fst_iterable(self.hor.get_mono_indexes()), 0)
