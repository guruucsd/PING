"""
Similarity matrix comparisons PING data.
"""
import copy
import sys
from functools import partial

from matplotlib import pyplot as plt

from ping.analysis.similarity import (compare_similarity_matrices,
                                      compute_similarity_matrices,
                                      visualize_similarity_matrices)
from ping.data import PINGData
from ping.utils import filter_dict
from research.asymmetry import is_ai_key
from research.data import get_all_data
from research.grouping import parse_filter_args


def do_usage(args):
    print("\nUsage for %s:" % args[0])
    print("\t%s [action] [measure]" % args[0])
    print("\n\taction: 'display' or 'launch'")
    print("\t\tdisplay: show results from a previous run")
    print("\t\tlaunch: launch a new GWAS run")
    print("\tmeasure: any measure from the PING database, or measure uploaded via upload.py")
    print("\t\tThe measure will be regressed against variation in genes at each SNP.")
    print("\t\tAge_At_ImgExam will be used as the covariate.")
    print("\nExamples:")
    print("\t%s launch MRI_cort_area_ctx_total_LH_PLUS_RH" % args[0])
    print("\t\tThis launches a GWAS to search for genetic variation as a function of total cortical area.")
    print("\t%s display MRI_cort_area_ctx_total_LH_PLUS_RH" % args[0])
    print("\t\tThis downloads and displays the top 200 SNPs related to total cortical area.")


if __name__ == '__main__':
    filter_args = parse_filter_args(sys.argv[1:])

    # Get data and compute matrices
    ai_group_fn = {
        'ai': lambda key: is_ai_key(key) and '_TOTAL' not in key}

    sim_dict = dict()
    for p in filter_args['prefixes']:
        prefix_filt_fn = partial(lambda k, v, p: k.startswith(p), p=p)

        p_data = get_all_data()
        p_data.filter(prefix_filt_fn)

        sim_dict[p] = compute_similarity_matrices(p_data.data_dict, filt_fns=ai_group_fn)

        # Compare and visualize
        compare_similarity_matrices(sim_dict[p])
        ax = visualize_similarity_matrices(sim_dict[p])
        ax.set_title('%s: %s' % (p, ax.get_title()))
    plt.show()
