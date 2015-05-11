"""
Similarity matrix comparisons PING data.
"""
import copy
from functools import partial

from matplotlib import pyplot as plt

from ping.analysis.similarity import (compare_similarity_matrices,
                                      compute_similarity_matrices,
                                      visualize_similarity_matrices)
from ping.data import PINGData
from ping.utils import filter_dict
from research.asymmetry import is_ai_key
from research.computed_measures import get_all_data


# Get data and compute matrices
prefixes = PINGData.IMAGING_PREFIX
ai_group_fn = {
    'ai': lambda key: is_ai_key(key) and '_TOTAL' not in key}

sim_dict = dict()
for p in prefixes:
    prefix_filt_fn = partial(lambda k, v, p: k.startswith(p), p=p)

    p_data = get_all_data()
    p_data.filter(prefix_filt_fn)

    sim_dict[p] = compute_similarity_matrices(p_data.data_dict, filt_fns=ai_group_fn)

    # Compare and visualize
    compare_similarity_matrices(sim_dict[p])
    ax = visualize_similarity_matrices(sim_dict[p])
    ax.set_title('%s: %s' % (p, ax.get_title()))
plt.show()
