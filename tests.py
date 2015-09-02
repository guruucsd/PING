from matplotlib import pyplot as plt

from similarity import do_similarity
from snps import do_snps


class TestWithNonblockingPlots(object):
    def setUp(self):
        plt.ion()


class TestSimilarity(TestWithNonblockingPlots):
    def test_similarity(self):
        do_similarity('MRI_cort_area.ctx', 'partial-correlation', 'Left Hemisphere')


class TestSnps(TestWithNonblockingPlots):
    # Only test anonymous functions
    def test_snp_view(self):
        do_snps('view', 'STK31')
