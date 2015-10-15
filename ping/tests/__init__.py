"""
"""
import os
import subprocess
import sys
import warnings

import matplotlib
from nose.tools import assert_in, assert_true, assert_raises, assert_equal
from six import StringIO, string_types


class TestWithNonblockingPlots(object):
    def setUp(self):
        matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!
        from matplotlib import pyplot as plt
        plt.ion()


class TestWithCaptureStdout(object):
    def setUp(self):
        # setup the environment
        self.backup, sys.stdout = sys.stdout, StringIO()

    def tearDown(self):
        sys.stdout.close()  # close the stream
        sys.stdout = self.backup  # restore original stdout

    def exec_pyfile(self, py_cmd):
        py_exec = sys.executable
        if isinstance(py_cmd, string_types):
            cmd = '%s %s' % (py_exec, py_cmd)
        elif isinstance(py_cmd, list):
            cmd = [py_exec] + py_cmd
        return self.exec_cmd(cmd)

    def exec_cmd(self, cmd):
        if isinstance(cmd, string_types):
            cmd = cmd.split(' ')

        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        return stdout.decode(), stderr.decode()


class TestWithGoodies(TestWithNonblockingPlots, TestWithCaptureStdout):
    def setUp(self):
        TestWithNonblockingPlots.setUp(self)
        TestWithCaptureStdout.setUp(self)


class TestExport(TestWithGoodies):
    def test_export_usage(self):
        warnings.warn('export is not using argparse.', DeprecationWarning)
        usage_text, stderr = self.exec_pyfile('export.py 1 2 3 4 5 6 7 8 9 10 11 12')
        assert_in("export.py", usage_text, 'usage with many parameters=%s' % usage_text)
        assert_equal(stderr, "", 'stderr should be empty; instead=="%s"' % stderr)

    def test_export(self):
        from export import EXPORTED_PING_SPREADSHEET
        if os.path.exists(EXPORTED_PING_SPREADSHEET):
            os.remove(EXPORTED_PING_SPREADSHEET)
        _, stderr = self.exec_pyfile('export.py')
        assert_true(os.path.exists(EXPORTED_PING_SPREADSHEET))
        assert_equal(stderr, "", 'stderr should be empty; instead=="%s"' % stderr)


class TestGrouping(TestWithGoodies):

    def test_grouping_usage(self):
        stdout, usage_text = self.exec_pyfile('grouping.py')
        assert_equal(stdout, "", 'stdout should be empty; instead=="%s"' % stdout)
        assert_in("grouping.py", usage_text, 'usage with no parameters')

        stdout, usage_text = self.exec_pyfile('grouping.py 1 2 3 4 5 6 7 8 9 10 11 12')
        assert_in("grouping.py", usage_text, 'usage with many parameters')
        assert_equal(stdout, "", 'stdout should be empty; instead=="%s"' % stdout)

    # def test_grouping(self):
    #     _, stderr = self.exec_pyfile('grouping.py MRI_cort_area.ctx Gender')
    #     assert_equal(stderr, "", 'Stderr should be empty; instead=="%s"' % stderr)


class TestGwas(TestWithGoodies):

    def test_gwas_usage(self):
        stdout, usage_text = self.exec_pyfile('gwas.py')
        assert_in("gwas.py", usage_text, 'usage with no parameters')
        assert_equal(stdout, "", 'stdout should be empty; instead=="%s"' % stdout)

        stdout, usage_text = self.exec_pyfile('gwas.py 1 2 3 4 5 6 7 8 9 10 11 12')
        assert_in("gwas.py", usage_text, 'usage with many parameters')
        assert_equal(stdout, "", 'stdout should be empty; instead=="%s"' % stdout)


class TestScatter(TestWithGoodies):

    def test_scatter_usage(self):
        stdout, usage_text = self.exec_pyfile('scatter.py')
        assert_in("scatter.py", usage_text, 'usage with no parameters')
        assert_equal(stdout, "", 'stdout should be empty; instead=="%s"' % stdout)

        stdout, usage_text = self.exec_pyfile('scatter.py 1 2 3 4 5 6 7 8 9 10 11 12')
        assert_in("scatter.py", usage_text, 'usage with many parameters')
        assert_equal(stdout, "", 'stdout should be empty; instead=="%s"' % stdout)

    # def test_scatter(self):
    #     _, stderr = self.exec_pyfile('scatter.py MRI_cort_area.ctx AI:mean AI:std LH_PLUS_RH:mean')
    #     assert_equal(stderr, "", 'stderr should be empty; instead=="%s"' % stderr)


class TestSimilarity(TestWithGoodies):

    def test_similarity_usage(self):
        stdout, usage_text = self.exec_pyfile('similarity.py')
        assert_in("similarity.py", usage_text, 'usage with no parameters')
        assert_equal(stdout, "", 'stdout should be empty; instead=="%s"' % stdout)

        stdout, usage_text = self.exec_pyfile('similarity.py 1 2 3 4 5 6 7 8 9 10 11 12')
        assert_in("similarity.py", usage_text, 'usage with many parameters')
        assert_equal(stdout, "", 'stdout should be empty; instead=="%s"' % stdout)

    # def test_similarity(self):
    #     _, stderr = self.exec_pyfile(['similarity.py', 'MRI_cort_area.ctx', 'partial-correlation', 'Asymmetry Index'])
    #     assert_equal(stderr, "", 'stderr should be empty; instead=="%s"' % stderr)


class TestSnps(TestWithGoodies):

    def test_snp_usage(self):
        stdout, usage_text = self.exec_pyfile('snps.py')
        assert_in("snps.py", usage_text, 'usage with no parameters')
        assert_equal(stdout, "", 'stdout should be empty; instead=="%s"' % stdout)

        stdout, usage_text = self.exec_pyfile('snps.py 1 2 3 4 5 6 7 8 9 10 11 12')
        assert_in("snps.py", usage_text, 'usage with many parameters')
        assert_equal(stdout, "", 'stdout should be empty; instead=="%s"' % stdout)

    def test_snp_view(self):
        _, stderr = self.exec_pyfile('snps.py view STK31')
        assert_equal(stderr, "", 'stderr should be empty; instead=="%s"' % stderr)
