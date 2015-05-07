"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import os

import csv
import matplotlib.pyplot as plt
import pandas
import numpy as np
import scipy

from export_measures import get_all_derived_data
from ping import col2prop, load_PING_data, get_twohemi_keys, get_asymmetry_index
from utils import do_and_plot_regression


def skip_key(key):
    return ('fuzzy' in key or
            key == 'MRI_cort_thick_ctx_mean_AI')


def skip_pairing(key1, key2):
    return (('AllFib' in key1 and 'AllFib' in key2) or
            ('AllFib' in key1 and 'DTI_' in key2) or
            ('AllFib' in key2 and 'DTI_' in key1) or
            np.any([t in key1 and t in key2 for t in ['SLF_', 'SCS_', '_Fx']]))


def find_one_relationship(all_data, key1, key2, covariates=[],
                          rsq_thresh=0., plot=False):
    print key1, key2, covariates

    # Limit to data without nan
    idx = np.ones(all_data[key1].shape, dtype=bool)
    for key in [key1, key2] + covariates:
        idx = np.logical_and(idx, np.logical_not(np.isnan(all_data[key])))
    if not np.any(idx):
        return None

    # Construct the data matrices
    data1 = all_data[key1][idx]
    data2 = all_data[key2][idx]
    X = [np.ones((idx.sum(),))]
    for key in [key2] + covariates:
        X.append(all_data[key][idx])

    # Do the regression, remove covariates.
    bestfit = np.linalg.lstsq(np.asarray(X).T, data1)[0]
    for ci, dat in enumerate(X[2:]):  # covariates
        data2 -= bestfit[ci + 2] * dat

    # Now, redo the regression on the residuals
    m, b, r, p, err = scipy.stats.linregress(data1, data2)
    assert np.logical_not(np.isnan(r)), 'WTF, nan?'
    if r**2 < rsq_thresh:
        return None

    key = '%s vs. %s' % (key1, key2)

    # Small plot
    if plot:
        xlims = np.asarray([np.min(data1), np.max(data1)])
        ylims = np.asarray([np.min(data2), np.max(data2)])
        xlims += (np.diff(xlims) / 10. * np.asarray([-1, 1]))  # add space
        ylims += (np.diff(ylims) / 10. * np.asarray([-1, 1]))  # add space

        # Grab / create the axis handle.
        if isinstance(plot, bool):
            ax = plt.figure().gca()
        else:
            ax = plot

        # Plot axes
        ax.hold(True)
        if xlims[0] < 0 and xlims[1] > 0:
            ax.plot([0, 0], ylims, 'k--')
        if ylims[0] < 0 and ylims[1] > 0:
            ax.plot(xlims, [0, 0], 'k--')
        ax.scatter(data1, data2)
        ax.plot(xlims, b + m * xlims, 'r', linewidth=3.0)
        x_mean = data1.mean()
        ax.plot(x_mean, b + m * x_mean, 'g*', markersize=20.0)

        # Metadata
        ax.set_title("Significant at %.2e (r=%.3f)" % (p, r))
        ax.set_xlabel(key1)
        ax.set_ylabel(key2)
        ax.set_xlim(tuple(xlims))
        ax.set_ylim(tuple(ylims))

    return key, p, r


def search_all_pairwise(all_data):
    """For each pair of variables, look for a significant regression slope."""
    results = []
    keys = list(set(all_data.keys()) - set(('SubjID',)))
    for ki in range(len(keys)):
        key1 = keys[ki]
        if skip_key(key1):
            continue

        for ii in range(ki + 1, len(keys)):
            key2 = keys[ii]
            if skip_key(key2) or skip_pairing(key1, key2):
                continue
            
            result = find_one_relationship(all_data, key1, key2, rsq_thresh=0.10)
            if result is not None:
                results.append(result)

    # Now, output the sorted result.
    for key, p, r in sorted(results, lambda v1, v2: int(10000 * (abs(v1[2]) - abs(v2[2])))):
        print "Significant at %.2e (r=%.3f): %s" % (p, r, key)



def search_all_vs_one(all_data, key, rsq_thresh=0., covariates=[], plot=False):
    """For each pair of variables, look for a significant regression slope."""

    results = []
    all_keys = list(set(all_data.keys()) - set(('SubjID',)))
    for all_key in all_keys:
        if not all_key.endswith('_AI'):
            continue
        try:
            result = find_one_relationship(all_data,
                                           key,
                                           all_key,
                                           covariates=covariates,
                                           rsq_thresh=rsq_thresh,
                                           plot=plot)
        except TypeError as te:
            # Happens when the data type is string or other non-regressable.
            continue
        if result is not None:
            results.append(result)

    # Now, output the sorted result.
    for all_key, p, r in sorted(results, lambda v1, v2: int(10000 * (v1[2] - v2[2]))):
        print "Significant at %.2e (r=%.3f): %s" % (p, r, all_key)



def print_legend():
    legend = {
        'CGH': "Parahippocampal Cingulum",
        'CST': "Corticospinal Tract",
        'IFO': "Inferior-Fronto-Occipital Fasiculus",
        "IFSFC": "Inferior Frontal Superior Frontal Cortex",
        "SCS": "Superior Corticostriate",
        "SLF": "Superior Longitudinal Fasiculus"}

    for key, val in legend.items():
        print '%5s: %s' % (key, val)

if __name__ == '__main__':
    all_data = get_all_derived_data(prefix=['MRI_cort_area', 'MRI_cort_thick',
                                            'MRI_subcort_vol', 'DTI_fiber_vol'])
    # search_all_pairwise(all_data)

    # search_all_vs_one(all_data, key='MRI_cort_area_ctx_total_LH_PLUS_RH', rsq_thresh=0.0015, plot=False)
    # search_all_vs_one(all_data, key='MRI_cort_area_TOTAL_AI', rsq_thresh=0.005, plot=True)
    search_all_vs_one(all_data, key='MRI_cort_area_ctx_entorhinal_AI', rsq_thresh=0.005)
    # search_all_vs_one(all_data, key='MRI_cort_area_ctx_frontalpole_AI',
    #                  covariates=['MRI_cort_area_ctx_total_LH_PLUS_RH', 'MRI_cort_thick_ctx_frontalpole_AI'],
    #                  rsq_thresh=0.005, plot=True)
    print_legend()
    plt.show()