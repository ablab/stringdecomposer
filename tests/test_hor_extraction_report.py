import os
import unittest

from sd.hor.hor_extraction_parser import HORExtractionReport

this_dirname = os.path.dirname(os.path.realpath(__file__))
test_data_dir = os.path.join(this_dirname, os.path.pardir, 'test_data')
monomers_fn = os.path.join(test_data_dir, 'DXZ1_star_monomers.fa')
sequences_fn = os.path.join(test_data_dir, 'read.fa')
sd_report_fn = os.path.join(test_data_dir, 'final_decomposition_fc89af8.tsv')
hor_report_fn = os.path.join(test_data_dir,
                             'cenX_hor_decomposition_f34264a.tsv')


class TestHORExtractionReport(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.hor_report = HORExtractionReport(hor_report_fn=hor_report_fn,
                                              sd_report_fn=sd_report_fn,
                                              monomers_fn=monomers_fn,
                                              sequences_fn=sequences_fn)
        self.horstring_set = self.hor_report.horstring_set
        self.s_id = "774d5aae-1eda-432b-a14f-3c18dad0b36b|runid=e9a32ee0ac1367543dc685eadc375ac1c17cb2b5|sampleid=CHM13_1|read=3871|ch=392|start_time=2018-05-01T21:31:27Z"
        super(TestHORExtractionReport, self).__init__(*args, **kwargs)

    def test_horstring_set(self):
        self.assertEqual(len(self.horstring_set), 1)
        horstring = self.horstring_set[self.s_id]
        self.assertEqual (horstring[0], (7,))
        self.assertEqual(horstring[10], tuple(range(12)))
        horstring.get_compact_form()
