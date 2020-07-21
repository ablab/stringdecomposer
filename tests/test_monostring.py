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


class TestMonostring(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        self.sd_report = SD_Report(sd_report_fn=sd_report_fn,
                                   monomers_fn=monomers_fn,
                                   sequences_fn=sequences_fn)
        self.monostring = fst_iterable(self.sd_report.monostring_set.values())
        self.seq_id = "774d5aae-1eda-432b-a14f-3c18dad0b36b|runid=e9a32ee0ac1367543dc685eadc375ac1c17cb2b5|sampleid=CHM13_1|read=3871|ch=392|start_time=2018-05-01T21:31:27Z"
        self.length = 557
        self.nucl_sequence = RC(read_bio_seq(sequences_fn))
        self.is_reversed = True
        self.corrections = {}
        self.fst_monoind = '?'
        self.snd_monoind = 6
        super(TestMonostring, self).__init__(*args, **kwargs)

    def test_properties(self):
        monostring = self.monostring
        self.assertEqual(monostring.seq_id, self.seq_id)
        self.assertEqual(len(monostring), self.length)
        self.assertEqual(len(monostring.raw_monostring), self.length)
        self.assertEqual(len(monostring.monoinstances), self.length)
        self.assertEqual(monostring.nucl_sequence, self.nucl_sequence)
        self.assertEqual(monostring.is_reversed, self.is_reversed)
        self.assertEqual(monostring.corrections, self.corrections)
        self.assertFalse(monostring.is_corrected())

    def test_getitem(self):
        self.assertEqual(self.monostring[0], self.fst_monoind)
        self.assertEqual(self.monostring[1], self.snd_monoind)

    def test_monomer_classification(self):
        classification = \
            self.monostring.classify_monomerinstances_by_monoindex()
        self.assertEqual(len(classification[0]), 46)
        self.assertEqual(len(self.monostring.
                             get_monomerinstances_by_monoindex(0)),
                         46)

    def test_get_identities(self):
        self.assertAlmostEqual(self.monostring.get_identities()[0], 0.6433)

    def test_get_perc_reliable(self):
        perc = 0.9964093
        self.assertAlmostEqual(self.monostring.get_perc_reliable(), perc)
        self.assertAlmostEqual(self.monostring.get_perc_unreliable(), 1-perc)

    def test_get_perc_strand(self):
        self.assertAlmostEqual(self.monostring.get_perc_forward_strand(), 1)
        self.assertAlmostEqual(self.monostring.get_perc_reverse_strand(), 0)

    def test_monoinstance(self):
        monoinstance = self.monostring.monoinstances[0]
        self.assertEqual(monoinstance.get_monoid(),
                         'F_5_DXZ1*_doubled/789_959/R')
        self.assertEqual(monoinstance.get_secmonoid(),
                         'K_10_DXZ1*_doubled/1639_1808/R')
        self.assertEqual(monoinstance.get_monoindex(), 5)
        self.assertEqual(monoinstance.get_secmonoindex(), 10)
        self.assertTrue(monoinstance.is_forward())
        self.assertFalse(monoinstance.is_reliable())
