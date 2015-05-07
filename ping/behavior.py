"""
"""

import numpy as np
import scipy
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

from ping import col2prop, get_asymmetry_index, load_PING_data, get_twohemi_keys
from similarity_matrices import get_all_data
from utils import do_and_plot_regression


def report_asymmetries():
    """Stale code moved from elsewhere."""
    # Report a summary of the asymmetries
    n_measures = len(good_keys)
    all_thresh = [1.0, 0.25, 0.10, 0.05, 0.01]
    for ti, thresh in enumerate(all_thresh[1:]):
        abs_mean = np.abs(mean_asymmetry)
        idx_above = np.logical_and(thresh <= abs_mean, abs_mean < all_thresh[ti])
        n_above = int(idx_above.astype(int).sum())
        print("%d (of %d) asymmetries > %.2f" % (n_above, n_measures, thresh))
        for hemi_sign, hemi_name in zip([1, -1], ['lh', 'rh']):
            hemi_idx = np.logical_and(np.sign(mean_asymmetry) == hemi_sign, idx_above)
            n_hemi = np.abs(hemi_idx.sum())
            print("\t%d (of %d) are %s" % (n_hemi, n_above, hemi_name))
            for idx in np.nonzero(hemi_idx.astype(int))[0]:
                print("\t\t%s (%.2f)" % (good_keys[idx], mean_asymmetry[idx]))


def show_asymmetry_by_component(prefix, grouping_prop_name='FDH_23_Handedness_Prtcpnt',
                                whiten=True, pc_threshhold=0.025):
    """ Grab all data matching prefix, split into groups by grouping_prop_name,
    then perform PCA and show the PCs."""

    data = get_all_data(prefix=prefix + ['TBX_', 'FDH_'])
    good_keys = get_twohemi_keys(prefix, data.keys())
    good_keys = np.asarray([key for key in good_keys if 'fuzzy' not in key])

    # Narrow to selected keys
    data_subset = []
    for key in good_keys:
        data_subset.append(get_asymmetry_index(data, key)) #/ data['MRI_cort_area_ctx_total_LH_PLUS_RH'])
    data_subset = np.asarray(data_subset)

    # Eliminate subjects with zero asymmetry in anything (likely artifacts)
    nan_idx = np.isnan(data_subset.sum(axis=0))
    asymm_mag = np.abs(data_subset).sum(axis=0)
    noasymm_idx = asymm_mag < asymm_mag.max() / 1000.
    good_subj_idx = np.logical_not(np.logical_or(nan_idx, noasymm_idx))
    data_subset = data_subset[:, good_subj_idx]
    print "Eliminated %d subjects with zero asymmetry, %d with a nan measure." % (
        sum(noasymm_idx), sum(nan_idx))
    # if whiten:
    # data_subset = scipy.stats.mstats.zscore(data_subset, axis=1)

    # Narrow to selected keys with mean(asymmetry_index) > 0.05
    mean_asymmetry = data_subset.std(axis=1)
    # asymmetric_idx = np.abs(mean_asymmetry) >= asymmetry_threshhold
    # data_subset = data_subset[asymmetric_idx]
    # good_keys = good_keys[asymmetric_idx]
    # mean_asymmetry = mean_asymmetry[asymmetric_idx]

    # Compute PCA and dump the results
    pca = PCA(whiten=whiten)
    pca.fit(data_subset.T)
    selected_idx = pca.explained_variance_ratio_ >= pc_threshhold
    for pi, pc in enumerate(pca.components_[selected_idx]):
        print "%2d: (%.2f):" % (pi, pca.explained_variance_ratio_[pi])
        sort_idx = np.argsort(np.abs(pc) * mean_asymmetry)
        for key, coeff, mag in zip(good_keys[sort_idx], pc[sort_idx], mean_asymmetry[sort_idx]):
            # if np.sign(mag) == hemi_sign:
            print "\t%s%.4f %-50s (%.3e / %.3e)" % (
                ' ' if mag*coeff >= 0 else '',
                mag*coeff,
                key,
                coeff,
                mag)

    # Now, take this factor and compare it to behavioral data
    # (or, should I throw the behavior in the PCA?)
    tbx_data = dict()
    for key in data.keys():
        if not key.startswith('TBX_'):
            continue
        if data[key].dtype.name in ['string', 'object']:
            continue
        tbx_data[key] = data[key][good_subj_idx]
    print 'Found %d TBX keys.' % len(tbx_data)

    for pi, pc in enumerate(pca.components_[selected_idx]):
        pc_projection = np.dot(pc, data_subset)
        for ti, (key, val) in enumerate(tbx_data.items()):
            cur_idx = np.logical_not(np.isnan(val))
            rval, pval = pearsonr(pc_projection[cur_idx], val[cur_idx])
            if pval < 0.05:
                print('%d: %-20s (n=%4d) r=%.4f, p=%.4f' % (pi, key, cur_idx.sum(), rval, pval))


    # Now, take this factor and compare it to hisory data
    # (or, should I throw the behavior in the PCA?)
    fdh_data = dict()
    for key in data.keys():
        if not key.startswith('FDH_'):
            continue
        if data[key].dtype.name in ['string', 'object']:
            continue
        fdh_data[key] = data[key][good_subj_idx]
    print 'Found %d FDH keys.' % len(fdh_data)

    for pi, pc in enumerate(pca.components_[selected_idx]):
        pc_projection = np.dot(pc, data_subset)
        for ti, (key, val) in enumerate(fdh_data.items()):
            cur_idx = np.logical_not(np.isnan(val))
            rval, pval = pearsonr(pc_projection[cur_idx], val[cur_idx])
            if pval < 0.05:
                print('%d: %-20s (n=%4d) r=%.4f, p=%.4f' % (pi, key, cur_idx.sum(), rval, pval))


if __name__ == '__main__':
    # 'MRI_cort_area', 'MRI_cort_thick', 'MRI_subcort_vol', 'DTI_fiber_vol'
    show_asymmetry_by_component(prefix=['DTI_fiber_vol'])
