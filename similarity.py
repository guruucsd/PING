"""
Similarity matrix comparisons PING data.
"""
import copy
import sys
from functools import partial

import numpy as np
from matplotlib import pyplot as plt

from ping.analysis.similarity import (compare_similarity_matrices,
                                      compute_similarity_matrices,
                                      visualize_similarity_matrices)
from ping.data import PINGData
from ping.utils import filter_dict
from research.asymmetry import is_ai_key
from research.data import get_all_data
from research.grouping import parse_filter_args


def do_usage(args, error_msg=None):
    if error_msg is not None:
        print("*** ERROR *** : %s" % error_msg)
    print("\nUsage: %s [prefix]" % args[0])
    print("\tCompare asymmetry correlation matrix with LH/RH structural covariance matrices.")
    print("\t\tprefix: [optional] comma-separated list of prefixes to include in the analysis.")


if __name__ != '__main__':
    pass

elif len(sys.argv) > 2:
    do_usage("Too many arguments.")

else:
    prefix = PINGData.IMAGING_PREFIX
    if len(sys.argv) == 2:
        prefix = sys.argv[1].split(',')
    prefix_filter_fn = lambda k, v: np.any([k.startswith(p) for p in prefix])

    # Get data and compute matrices
    ai_group_fn = {
        'ai': lambda key: is_ai_key(key) and '_TOTAL' not in key}

    p_data = get_all_data()
    p_data.filter(prefix_filter_fn)

    sim_dict, good_keys = compute_similarity_matrices(p_data.data_dict,
                                                      filt_fns=ai_group_fn)
    for k in good_keys:
        print(k)

    # Compare and visualize
    compare_similarity_matrices(sim_dict)
    ax = visualize_similarity_matrices(sim_dict, labels=good_keys)
    ax.get_figure().suptitle(', '.join(prefix))

    plt.show()
