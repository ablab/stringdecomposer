import os
import unittest

from sd.sd_parser.sd_parser import SD_Report


this_dirname = os.path.dirname(os.path.realpath(__file__))
test_data_dir = os.path.join(this_dirname, os.path.pardir, 'test_data')
monomers_fn = os.path.join(test_data_dir, 'DXZ1_star_monomers.fa')
sequences_fn = os.path.join(test_data_dir, 'read.fa')
sd_report_fn = os.path.join(test_data_dir, 'final_decomposition_fc89af8.tsv')


class TestSDReport(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.sd_report = SD_Report(sd_report_fn=sd_report_fn,
                                   monomers_fn=monomers_fn,
                                   sequences_fn=sequences_fn)
        self.db_size = 12
        super(TestSDReport, self).__init__(*args, **kwargs)

    def test_sd_report(self):
        self.assertEqual(self.sd_report.monomer_db.get_size(), self.db_size)
