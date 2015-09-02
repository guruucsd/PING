"""
"""
import os

from matplotlib import pyplot as plt
from nose.tools import assert_true

from export import do_export, EXPORTED_PING_SPREADSHEET
from grouping import do_grouping
from scatter import do_scatter
from similarity import do_similarity
from snps import do_snps


class TestWithNonblockingPlots(object):
    def setUp(self):
        plt.ion()


class TestExport(TestWithNonblockingPlots):
    def test_export(self):
        do_export()
        assert_true(os.path.exists(EXPORTED_PING_SPREADSHEET))


class TestGrouping(TestWithNonblockingPlots):
    def test_grouping(self):
        do_grouping('MRI_cort_area.ctx', 'Gender')


class TestScatter(TestWithNonblockingPlots):
    def test_scatter(self):
        do_scatter('MRI_cort_area.ctx', 'AI:mean', 'AI:std', 'LH_PLUS_RH:mean')


class TestSimilarity(TestWithNonblockingPlots):
    def test_similarity(self):
        do_similarity('MRI_cort_area.ctx', 'partial-correlation', 'Left Hemisphere')


class TestSnps(TestWithNonblockingPlots):
    # Only test anonymous functions
    def test_snp_view(self):
        do_snps('view', 'STK31')
