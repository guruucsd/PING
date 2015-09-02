from matplotlib import pyplot as plt

from scatter import do_scatter
from similarity import do_similarity
from snps import do_snps


class TestWithNonblockingPlots(object):
    def setUp(self):
        plt.ion()


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
