"""
"""
import sys

from nose.tools import assert_raises, assert_raises_regexp


def test_export():
    import export


def test_grouping():
    import grouping


def test_gwas():
    sys.argv = [sys.argv[0], 'display', 'MRI_cort_area_ctx_rh_frontalpole']
    try:
        import gwas
    except Exception as e:
        assert_true('id not found' in str(e))


def test_regress():
    # with assert_raises(NotImplementedError):
    #     import regress
    pass


def test_similarity():
    import similarity


def test_upload():
    sys.argv = sys.argv[0:1]
    import upload
