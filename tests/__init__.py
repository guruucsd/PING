"""
"""
import os
import sys

from matplotlib import pyplot as plt
from nose.tools import assert_in, assert_true
from six import StringIO

from export import do_export, EXPORTED_PING_SPREADSHEET
from gwas import do_gwas
from grouping import do_grouping
from scatter import do_scatter
from similarity import do_similarity
from snps import do_snps


class TestWithNonblockingPlots(object):
    def setUp(self):
        plt.ion()


class TestWithCaptureStdout(object):
    def setUp(self):
        # setup the environment
        self.backup, sys.stdout = sys.stdout, StringIO()

    def tearDown(self):
        sys.stdout.close()  # close the stream
        sys.stdout = self.backup  # restore original stdout


class TestWithGoodies(TestWithNonblockingPlots, TestWithCaptureStdout):
    def setUp(self):
        TestWithNonblockingPlots.setUp(self)
        TestWithCaptureStdout.setUp(self)


class TestExport(TestWithGoodies):

    def test_export_usage(self):
        do_export(*range(12))
        usage_text = sys.stdout.getvalue()  # release output
        assert_in("export.py", usage_text, usage_text)

    def test_export(self):
        do_export()
        assert_true(os.path.exists(EXPORTED_PING_SPREADSHEET))


class TestGrouping(TestWithGoodies):

    def test_grouping_usage(self):
        do_grouping(*range(12))
        usage_text = sys.stdout.getvalue()  # release output
        assert_in("grouping.py", usage_text, usage_text)

    def test_grouping(self):
        do_grouping('MRI_cort_area.ctx', 'Gender')


class TestGwas(TestWithGoodies):

    def test_gwas_usage(self):
        do_gwas(*range(12))
        usage_text = sys.stdout.getvalue()  # release output
        assert_in("gwas.py", usage_text, usage_text)

    def test_grouping(self):
        do_grouping('MRI_cort_area.ctx', 'Gender')


class TestScatter(TestWithGoodies):

    def test_export_usage(self):
        do_scatter(*range(12))
        usage_text = sys.stdout.getvalue()  # release output
        assert_in("scatter.py", usage_text, usage_text)

    def test_scatter(self):
        do_scatter('MRI_cort_area.ctx', 'AI:mean', 'AI:std', 'LH_PLUS_RH:mean')


class TestSimilarity(TestWithGoodies):

    def test_export_usage(self):
        do_similarity(*range(12))
        usage_text = sys.stdout.getvalue()  # release output
        assert_in("similarity.py", usage_text, usage_text)

    def test_similarity(self):
        do_similarity('MRI_cort_area.ctx', 'partial-correlation', 'Left Hemisphere')


class TestSnps(TestWithGoodies):

    def test_snp_usage(self):
        do_snps(*range(12))
        usage_text = sys.stdout.getvalue()  # release output
        assert_in("snps.py", usage_text, usage_text)

    # Only test anonymous functions
    def test_snp_view(self):
        do_snps('view', 'STK31')
