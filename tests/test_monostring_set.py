import os
import unittest

from sd.sd_parser.sd_parser import SD_Report
from sd.utils.bio import read_bio_seq, RC
from sd.utils.various import fst_iterable

this_dirname = os.path.dirname(os.path.realpath(__file__))
test_data_dir = os.path.join(this_dirname, os.path.pardir, 'test_data')
monomers_fn = os.path.join(test_data_dir, 'DXZ1_star_monomers.fa')
sequences_fn = os.path.join(test_data_dir, 'read.fa')
sd_report_fn = os.path.join(test_data_dir, 'final_decomposition_fc89af8.tsv')


class TestMonostringSet(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.sd_report = SD_Report(sd_report_fn=sd_report_fn,
                                   monomers_fn=monomers_fn,
                                   sequences_fn=sequences_fn)
        self.monostring_set = self.sd_report.monostring_set
        super(TestMonostringSet, self).__init__(*args, **kwargs)

    def test_len(self):
        self.assertEqual(len(self.monostring_set), 1)

    def test_monomer_classification(self):
        classification = \
            self.monostring_set.classify_monomerinstances_by_monoindex()
        self.assertEqual(len(classification[0]), 46)
        self.assertEqual(len(self.monostring_set.
                             get_monomerinstances_by_monoindex(0)),
                         46)
