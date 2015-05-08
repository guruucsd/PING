"""
Multivariate analyses of PING data
"""
import numpy as np
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

from .access import (load_PING_data, get_twohemi_keys,
                     get_fdh_data, get_tbx_data)
from .asymmetry import get_asymmetry_index, is_ai_prop_name
from .utils import do_and_plot_regression


class AsymmetryPCA(object):
    def __init__(self, whiten=True, pc_threshhold=0.05):
        self.pca = None
        self.data = None
        self.data_mat = None
        self.good_keys = None

        self.whiten = whiten
        self.pc_threshhold = pc_threshhold

    def fit(self, data):
        """Perform PCA on asymmetry indices.
        data is a data dictionary."""

        data_mat = []
        keys = [key for key in data.keys()
                if 'fuzzy' not in key and '_TOTAL_AI' not in key]
        good_keys = get_twohemi_keys(keys)
        if len(good_keys) == 0:
            # We got asymmetry indices directly
            good_keys = [key for key in keys if is_ai_prop_name(key)]
            for key in good_keys:
                data_mat.append(data[key])
        else:
            # We got raw data
            good_keys = np.asarray(keys)
            for key in good_keys:
                data_mat.append(get_asymmetry_index(data, key))

        good_keys = np.asarray(good_keys)
        data_mat = np.asarray(data_mat)
        print data_mat.shape

        # Eliminate subjects with zero asymmetry in anything (likely artifacts)
        key_nan_idx = np.isnan(data_mat).sum(1) >= 0.95 * data_mat.shape[1]
        good_key_idx = np.logical_not(key_nan_idx)
        good_keys = good_keys[good_key_idx]
        data_mat = data_mat[good_key_idx]

        asymm_mag = np.abs(data_mat).sum(0)
        noasymm_idx = asymm_mag < asymm_mag.max() / 1000.
        subj_nan_idx = np.isnan(data_mat.sum(0))  # Eliminate subjects
        good_subj_idx = np.logical_not(np.logical_or(subj_nan_idx,
                                                     noasymm_idx))
        data_mat = data_mat[:, good_subj_idx]

        print "Eliminated %d keys, %d subjects." % (
            sum(np.logical_not(good_key_idx)),
            sum(np.logical_not(good_subj_idx)))

        # Compute PCA and dump the results
        self.subj_ids = data['SubjID'][good_subj_idx]
        self.data = data
        self.data_mat = data_mat
        self.good_keys = good_keys

        self.pca = PCA(whiten=self.whiten)
        self.pca.fit(data_mat.T)

    def get_components(self):
        selected_idx = self.pca.explained_variance_ratio_ >= self.pc_threshhold
        return self.pca.components_[selected_idx]

    def get_projections(self):
        return np.dot(self.get_components(), self.data_mat)

    def report_asymmetry_loadings(self):
        mean_asymmetry = self.data_mat.mean(axis=1)

        for pi, pc in enumerate(self.get_components()):
            print "%2d: (%.2f):" % (pi, self.pca.explained_variance_ratio_[pi])
            sort_idx = np.argsort(np.abs(pc * mean_asymmetry))
            for key, coeff, mag in zip(self.good_keys[sort_idx], pc[sort_idx], mean_asymmetry[sort_idx]):
                # if np.sign(mag) == hemi_sign:
                print "\t%s%.4f %-50s (%.3e / %.3e)" % (
                    ' ' if mag*coeff >= 0 else '',
                    mag*coeff,
                    key,
                    coeff,
                    mag)

        return self.pca

    def report_behavior_correlations(self):
        # Now, take this factor and compare it to behavioral data
        # (or, should I throw the behavior in the PCA?)
        tbx_data = get_tbx_data(dict(zip(self.good_keys, self.data_mat)))
        print 'Found %d TBX keys.' % len(tbx_data)

        for pi, proj in enumerate(self.get_projections()):
            for ti, (key, val) in enumerate(tbx_data.items()):
                cur_idx = np.logical_not(np.isnan(val))
                rval, pval = pearsonr(proj[cur_idx], val[cur_idx])
                if pval < 0.05:
                    print('%d: %-20s (n=%4d) r=%.4f, p=%.4f' % (pi, key, cur_idx.sum(), rval, pval))

    def report_background_correlations(self):
        # Now, take this factor and compare it to hisory data
        # (or, should I throw the behavior in the PCA?)
        fdh_data = get_fdh_data(dict(zip(self.good_keys, self.data_mat)))
        print 'Found %d FDH keys.' % len(fdh_data)

        for pi, pc in enumerate(self.get_projections()):
            for ti, (key, val) in enumerate(fdh_data.items()):
                cur_idx = np.logical_not(np.isnan(val))
                rval, pval = pearsonr(proj[cur_idx], val[cur_idx])
                if pval < 0.05:
                    print('%d: %-20s (n=%4d) r=%.4f, p=%.4f' % (pi, key, cur_idx.sum(), rval, pval))
