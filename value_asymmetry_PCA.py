"""
"""

import numpy as np
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

from ping import col2prop, get_asymmetry_index, load_PING_data, get_twohemi_keys
from utils import do_and_plot_regression


def show_asymmetry_by_component(prefix, grouping_prop_name='FDH_23_Handedness_Prtcpnt', whiten=True, asymmetry_threshhold=0.05):
    """ Grab all data matching prefix, split into groups by grouping_prop_name,
    then perform PCA and show the PCs."""

    data = load_PING_data()
    good_keys = get_twohemi_keys(prefix, data.keys())

    # Narrow to selected keys
    data_subset = []
    for key in good_keys:
        data_subset.append(get_asymmetry_index(data, key))
    data_subset = np.asarray(data_subset)

    # Eliminate subjects with zero asymmetry in anything (likely artifacts)
    noasymm_idx = np.abs(data_subset).sum(axis=0) != 0.
    print("Removing %d subjects with ZERO asymmetry "
          "(over %d measures!)" % (int(np.logical_not(noasymm_idx).sum()),
                                   len(good_keys)))
    data_subset = data_subset[:, noasymm_idx]

    # Narrow to selected keys with mean(asymmetry_index) > 0.05
    mean_asymmetry = data_subset.mean(axis=1)
    asymmetric_idx = np.abs(mean_asymmetry) >= asymmetry_threshhold
    data_subset = data_subset[asymmetric_idx]
    good_keys = good_keys[asymmetric_idx]
    mean_asymmetry = mean_asymmetry[asymmetric_idx]

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

    # Compute PCA and dump the results
    pca = PCA(whiten=whiten)
    pca.fit(data_subset.T)
    selected_idx = pca.explained_variance_ratio_ >= 0.05
    for pi, pc in enumerate(pca.components_[selected_idx]):
        print "%2d: (%.2f):" % (pi, pca.explained_variance_ratio_[pi])
        for hemi_sign, hemi_name in zip([1, -1], ['lh', 'rh']):
            for key, coeff, mag in zip(good_keys, pc, mean_asymmetry):
                if np.sign(mag) == hemi_sign:
                    print "\t%s %s%.4f %s (%.2f)" % (hemi_name, ' ' if coeff >= 0 else '', coeff, key, 1000*mag*coeff)

    # Now, take this factor and compare it to behavioral data
    # (or, should I throw the behavior in the PCA?)
    tbx_props = filter(lambda k: 'TBX_' in k, data.keys())
    tbx_data = np.asarray(data[tbx_props].to_records())[noasymm_idx]
    for pi, pc in enumerate(pca.components_[selected_idx]):
        pc_projection = np.dot(pc, data_subset)
        for ti, tbx_prop in enumerate(tbx_props):
            row_data = np.asarray([r[ti] for r in tbx_data])
            if 'string' in row_data.dtype.name:
                continue
            good_subj = np.logical_not(np.isnan(row_data))
            rval, pval = pearsonr(pc_projection[good_subj], row_data[good_subj])
            if pval < 0.05:
                print('%d: %s %.4f, %.4f' % (pi, tbx_prop, rval, pval))


if __name__ == '__main__':
    show_asymmetry_by_component(prefix=['MRI_subcort_vol', 'MRI_cort_thick',
                                        'MRI_cort_area', 'DTI_fiber_vol'])
