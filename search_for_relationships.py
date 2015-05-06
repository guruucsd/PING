"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import os

import csv
import matplotlib.pyplot as plt
import numpy as np
import scipy
import pandas

from export_measures import compute_all_asymmetries
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


def find_one_relationship(all_data, key1, key2, rsq_thresh=0., plot=False):
    idx1 = np.logical_not(np.isnan(all_data[key1]))
    idx2 = np.logical_not(np.isnan(all_data[key2]))
    idx = np.logical_and(idx1, idx2)

    if not np.any(idx):
        return None

    m, b, r, p, err = scipy.stats.linregress(all_data[key1][idx],
                                             all_data[key2][idx])
    assert np.logical_not(np.isnan(r)), 'WTF, nan?'
    if r**2 < rsq_thresh:
        return None

    key = '%s vs. %s' % (key1, key2)

    # Small plot
    if plot:
        xlims = np.asarray([np.min(all_data[key1][idx]), np.max(all_data[key1][idx])])
        ylims = np.asarray([np.min(all_data[key2][idx]), np.max(all_data[key2][idx])])
        xlims += (np.diff(xlims) / 10. * np.asarray([-1, 1]))  # add space
        ylims += (np.diff(ylims) / 10. * np.asarray([-1, 1]))  # add space

        # Grab / create the axis handle.
        if isinstance(plot, bool):
            ax = plt.figure()
        else:
            ax = plot

        # Plot axes
        ax.hold(True)
        if xlims[0] < 0 and xlims[1] > 0:
            ax.plot([0, 0], ylims, 'k--')
        if ylims[0] < 0 and ylims[1] > 0:
            ax.plot(xlims, [0, 0], 'k--')
        ax.scatter(all_data[key1][idx], all_data[key2][idx])
        ax.plot(xlims, b + m * xlims, 'r', linewidth=2.0)
        ax.set_title("Significant at %.2e (r=%.3f)" % (p, r))
        ax.set_xlabel(key1)
        ax.set_ylabel(key2)
        ax.set_xlim(tuple(xlims))
        ax.set_ylim(tuple(ylims))

    return key, p, r


def find_all_relationships(all_data):
    """For each pair of variables, look for a significant regression slope."""
    def skip_key(key):
        return ('fuzzy' in key or
                key == 'MRI_cort_thick_ctx_mean_AI')

    def skip_pairing(key1, key2):
        return (('AllFib' in key1 and 'AllFib' in key2) or
                ('AllFib' in key1 and 'DTI_' in key2) or
                ('AllFib' in key2 and 'DTI_' in key1) or
                np.any([t in key1 and t in key2 for t in ['SLF_', 'SCS_', '_Fx']]))

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
    all_data = compute_all_asymmetries(prefix=['MRI_cort_area', 'MRI_cort_thick',
                                               'MRI_subcort_vol', 'DTI_fiber_vol'])
    find_all_relationships(all_data)
    print_legend()
    plt.show()