"""
Build similarity matrices for cortical area (left, right) and asymmetry
"""
from collections import OrderedDict

import numpy as np
import scipy.spatial
import scipy.stats
from matplotlib import pyplot as plt

from .partial_corr import partial_corr
from .stats import pdist_nan
from ..utils.plotting import plot_symmetric_matrix_as_triangle


def is_bad_key(key):
    if np.any([substr in key.lower().replace('_', '.')
               for substr in ['.vent', 'fuzzy', 'bankssts', '.total']]):
        return True
    elif np.any([substr in key.lower().replace('_', '.')
                for substr in ['.mean', '.white.matter', '.cortex']]):
        return True
    elif np.any([substr in key.lower().replace('_', '.')
                for substr in ['allfib', '_slf', '_scs', '_fxcut']]):
        return True
    return False


def get_good_keys(all_data, filter_fn, sort_fn=sorted):
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
        elif all_data[key][~np.isnan(all_data[key])].std() == 0:
            rej_reason.append('std=0')
            continue  # Data without variation
        elif is_bad_key(key):
            rej_reason.append('substr')
            continue
        good_keys.append(key)

    if len(good_keys) == 0:
        raise ValueError("Filtered out all fields!")

    if sort_fn is not None:
        return sort_fn(good_keys)
    else:
        return good_keys


def build_similarity_vector(all_data, good_keys=None, filter_fn=None,
                            standardize=False, sort_fn=sorted,
                            metric='correlation'):
    """
    """
    # Filter keys
    if not good_keys:
        good_keys = sort_fn(get_good_keys(all_data, filter_fn))

    # Convert data dictionary into a data matrix
    data_mat = []
    for key in good_keys:
        data_mat.append(all_data[key])
    data_mat = np.asarray(data_mat)

    if metric != 'partial-correlation':
        dist_mat = pdist_nan(data_mat, metric=metric, standardize=standardize)
        corr_mat = 1 - dist_mat

    else:
        # Remove subjects with any nan, and standardize values
        # bad_idx = np.isnan(data_mat.sum(axis=0))
        # good_idx = np.logical_not(bad_idx)
        # assert good_idx.sum() > 0, "You rejected all data!"

        # data_mat = data_mat[:, good_idx]
        if standardize:
            for ri in range(data_mat.shape[0]):
                data_mat[ri] = scipy.stats.mstats.zscore(
                    data_mat[~np.isnan(data_mat[ri])])
            # assert np.all(np.abs(data_mat.sum(axis=1)) < 1E-4)
        # print("Found %d keys; removed %d subjects w/ missing data." % (
            # data_mat.shape[1], bad_idx.sum()))
        corr_mat = partial_corr(data_mat.T, verbose=1)
        corr_mat[np.isnan(corr_mat)] = 0  # some values with std=0
    assert not np.isnan(corr_mat.sum())

    return corr_mat, good_keys


def compare_similarity_matrices(vec1, vec2, metric='correlation'):

    # Make sure both similarity matrices are in vector form.
    if metric == 'correlation':
        # Convert to vector form
        if vec1.ndim == 2:
            vec1[np.eye(vec1.shape[0], dtype=bool)] = 0.
            vec1 = scipy.spatial.distance.squareform(vec1, 'tovector')
        if vec2.ndim == 2:
            vec2[np.eye(vec2.shape[0], dtype=bool)] = 0.
            vec2 = scipy.spatial.distance.squareform(vec2, 'tovector')

        # Now compare.
        return scipy.stats.pearsonr(vec1, vec2)

    elif metric == 'norm-ish':
        # Convert to matrix form
        if vec1.ndim == 1:
            vec1 = scipy.spatial.distance.squareform(vec1)
        if vec2.ndim == 1:
            vec2 = scipy.spatial.distance.squareform(vec2)
        return (np.trace(np.dot(vec1 - vec2, vec1 - vec2)), np.nan)


def compute_similarity_vectors(data, filt_fns=None, **kwargs):
    if filt_fns is None:
        # Default filter: one group of everything!
        filt_fns = dict(all=lambda k: True)

    sim_dict = OrderedDict()
    sim_keys = OrderedDict()
    for mat_type, filt_fn in filt_fns.items():
        print("Computing similarity matrix for %s" % (mat_type))

        sim_mat, good_keys = build_similarity_vector(data,
                                                     filter_fn=filt_fn,
                                                     **kwargs)
        sim_dict[mat_type] = sim_mat
        sim_keys[mat_type] = good_keys

    assert len(np.unique([len(vals) for vals in sim_dict.values()])) == 1

    return sim_dict, sim_keys


def compare_similarity_vectors(sim_dict):
    # 2. Compare similarity matrices.
    compare_keys = list(sim_dict.keys())
    n_keys = len(compare_keys)

    mat_compare_mat = np.zeros((n_keys * (n_keys - 1) / 2,))
    mat_idx = 0
    for ki, key1 in enumerate(compare_keys):
        for kj in range(ki + 1, n_keys):
            key2 = compare_keys[kj]
            r, pval = compare_similarity_matrices(sim_dict[key1],
                                                  sim_dict[key2])
            print("%s vs. %s: r**2=%.3f (p=%.3f)" % (
                key1, key2, r**2, pval))
            mat_compare_mat[mat_idx] = r
            mat_idx += 1


def visualize_similarity_matrices(sim_dict, labels=None, class_labels=None, dynamic_color=True):
    # Visualize similarity matrices
    compare_keys = list(sim_dict.keys())
    n_keys = len(compare_keys)
    fh = plt.figure(figsize=(18, 10))

    for ki, key in enumerate(compare_keys):
        vmin, vmax = -1, 1
        mat = sim_dict[key]
        if mat.ndim == 1:
            mat = scipy.spatial.distance.squareform(mat)
        if dynamic_color:
            mat[np.eye(mat.shape[0], dtype=bool)] = 0.
            vval = np.max(np.abs([mat.min(), mat.max()]))
            vmin, vmax = np.asarray([-1, 1]) * vval

        ax = fh.add_subplot(1, n_keys, ki + 1)
        plot_symmetric_matrix_as_triangle(mat, ax=ax,
                                          vmin=vmin, vmax=vmax,
                                          labels=labels if ki == 0 else None,
                                          class_labels=class_labels if ki == 0 else None)
        plt.tight_layout(h_pad=5)
        ax.set_title(key, fontsize=18)

    return ax
