"""
Similarity matrix comparisons PING data.
"""
from matplotlib import pyplot as plt

from data import get_all_data
from ping.similarity import (compare_similarity_matrices,
                             compute_similarity_matrices,
                             visualize_similarity_matrices)
from ping.utils import filter_dict


# Get data and compute matrices
prefixes = ['MRI_cort_area', 'MRI_cort_thick', 'MRI_subcort_vol', 'DTI_fiber_vol']

sim_dict = dict()
for p in prefixes:
    p_data = get_all_data(prefix=[p])
    p_data = filter_dict(p_data, lambda key, val: key.startswith(p))
    sim_dict[p] = compute_similarity_matrices(p_data)

    # Compare and visualize
    compare_similarity_matrices(sim_dict[p])
    ax = visualize_similarity_matrices(sim_dict[p])
    ax.set_title('%s: %s' % (p, ax.get_title()))
plt.show()
