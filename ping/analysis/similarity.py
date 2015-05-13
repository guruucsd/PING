"""
Build similarity matrices for cortical area (left, right) and asymmetry
"""
import copy
import numpy as np
import scipy.spatial
import scipy.stats
from matplotlib import pyplot as plt

from .partial_corr import partial_corr
from ..data import which_hemi
from ..utils.plotting import plot_symmetric_matrix_as_triangle


def get_good_keys(all_data, filter_fn):
    good_keys = []
    rej_reason = []
    for key in all_data.keys():
        if not filter_fn(key):
            rej_reason.append('filter')
            continue  # Bad key
        elif all_data[key].dtype.name in ['string', 'object']:
            rej_reason.append('dtype')
            continue
        elif np.isnan(all_data[key]).sum() == len(all_data[key]):
            rej_reason.append('all nan')
            continue  # All data is nan!
        elif all_data[key][np.logical_not(np.isnan(all_data[key]))].std() == 0:
            rej_reason.append('std=0')
            continue  # Data without variation
        elif np.any([substr in key.lower()
                    for substr in ['.vent', 'fuzzy', 'total', '.mean', '.white.matter', '.cortex', 'allfib']]):
            rej_reason.append('substr')
            continue  # Remove ventricles
        good_keys.append(key)
    assert len(good_keys) > 0
    return sorted(good_keys)


def build_similarity_matrix(all_data, good_keys=None, filter_fn=None, standardize=False):
    """
    """
    # Filter keys
    if not good_keys:
        good_keys = sorted(get_good_keys(all_data, filter_fn))

    # Convert data dictionary into a data matrix
    data_mat = []
    for key in good_keys:
        data_mat.append(all_data[key])
    data_mat = np.asarray(data_mat)

    # Remove subjects with any nan, and standardize values
    bad_idx = np.isnan(data_mat.sum(axis=0))
    good_idx = np.logical_not(bad_idx)
    data_mat = data_mat[:, good_idx]
    if standardize:
        data_mat = scipy.stats.mstats.zscore(data_mat, axis=1)
        assert np.all(np.abs(data_mat.sum(axis=1)) < 1E-4)
    print("Found %d keys; removed %d subjects w/ missing data." % (
        data_mat.shape[1], bad_idx.sum()))

    # Compute a correlation matrix
    dist_mat = scipy.spatial.distance.pdist(data_mat, 'correlation')
    corr_mat = partial_corr(data_mat.T)

    assert not np.isnan(corr_mat.sum())

    return corr_mat, good_keys


def compare_similarity_vectors(vec1, vec2, metric='correlation'):

    # Make sure both similarity matrices are in vector form.
    if metric == 'correlation':
        if len(vec1.shape) == 2:
            vec1[np.eye(vec1.shape[0], dtype=bool)] = 0.
            vec1 = scipy.spatial.distance.squareform(vec1, 'tovector')
        if len(vec2.shape) == 2:
            vec2[np.eye(vec2.shape[0], dtype=bool)] = 0.
            vec2 = scipy.spatial.distance.squareform(vec2, 'tovector')

        # Now compare.
        return scipy.stats.pearsonr(vec1, vec2)

    elif metric == 'norm-ish':
        if len(vec1.shape) == 1:
            vec1 = scipy.spatial.distance.squareform(vec1)
        if len(vec2.shape) == 1:
            vec2 = scipy.spatial.distance.squareform(vec2)
        return (np.trace(np.dot(vec1 - vec2, vec1 - vec2)), np.nan)


def compute_similarity_matrices(data, filt_fns=None):

    # Add default filters
    all_filters = {
        'left': lambda key: which_hemi(key) == 'lh',
        'right': lambda key: which_hemi(key) == 'rh'}
    if filt_fns is not None:
        all_filters.update(filt_fns)

    sim_dict = dict()

    # 1. Compute similarity matrices
    for mat_type, filt_fn in all_filters.items():
        print("Computing similarity matrix for %s" % (mat_type))

        sim_mat, good_keys = build_similarity_matrix(data,
                                                     filter_fn=filt_fn,
                                                     standardize=False)
        sim_dict[mat_type] = sim_mat
        print([v.shape for v in sim_dict.values()])
    return sim_dict, good_keys


def compare_similarity_matrices(sim_dict):
    # 2. Compare similarity matrices.
    compare_keys = list(sim_dict.keys())
    n_keys = len(compare_keys)

    mat_compare_mat = np.zeros((n_keys * (n_keys - 1) / 2,))
    mat_idx = 0
    for ki, key1 in enumerate(compare_keys):
        for kj in range(ki + 1, n_keys):
            key2 = compare_keys[kj]
            r, pval = compare_similarity_vectors(sim_dict[key1],
                                                 sim_dict[key2])
            print("%s vs. %s: r**2=%.3f (p=%.3f)" % (
                key1, key2, r**2, pval))
            mat_compare_mat[mat_idx] = r
            mat_idx += 1


def visualize_similarity_matrices(sim_dict, labels=None, dynamic_color=True):
    # Visualize similarity matrices
    compare_keys = list(sim_dict.keys())
    n_keys = len(compare_keys)

    fh = plt.figure(figsize=(16, 6))
    for ki, key in enumerate(compare_keys):
        vmin, vmax = -1, 1
        if dynamic_color:
            vval = np.max(np.abs([sim_dict[key].min(), sim_dict[key].max()]))
            vmin, vmax = np.asarray([-1, 1]) * vval

        ax = fh.add_subplot(1, n_keys, ki + 1)
        plot_symmetric_matrix_as_triangle(sim_dict[key], ax=ax,
                                          vmin=vmin, vmax=vmax,
                                          lbls=labels if ki == 0 else None)
        ax.set_title(key)

    return ax
