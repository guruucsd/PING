"""
Multivariate analyses of PING data
"""
import copy

import numpy as np
from scipy.stats import pearsonr
from sklearn.decomposition import PCA

from .asymmetry import get_asymmetry_index, is_ai_key
from ping.utils import do_and_plot_regression


class AsymmetryPCA(object):
    def __init__(self, whiten=True, pc_threshhold=0.05):
        self.pca = None
        self.data = None
        self.data_mat = None
        self.good_keys = None

        self.whiten = whiten
        self.pc_threshhold = pc_threshhold

    def fit(self, data, verbose=1):
        """Perform PCA on asymmetry indices.
        data is a data dictionary."""

        data = copy.deepcopy(data)
        data.filter(lambda k, v: 'fuzzy' not in k)
        data.filter(lambda k, v: '_TOTAL_AI' not in k)
        self.data = data
        good_keys = self.data.get_twohemi_keys()

        data_mat = []
        if len(good_keys) == 0:
            # We got asymmetry indices directly
            good_keys = [key for key in data.data_dict.keys() if is_ai_key(key)]
            for key in good_keys:
                data_mat.append(data.data_dict[key])
        else:
            # We got raw data
            good_keys = np.asarray(list(data.data_dict.keys()))
            for key in good_keys:
                data_mat.append(get_asymmetry_index(data.data_dict, key))

        good_keys = np.asarray(good_keys)
        data_mat = np.asarray(data_mat)

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

        if verbose >= 1:
            print("Eliminated %d keys, %d subjects." % (
                sum(np.logical_not(good_key_idx)),
                sum(np.logical_not(good_subj_idx))))

        print('PCA data matrix size: %s' % str(data_mat.shape))

        # Compute PCA and dump the results
        self.subj_ids = data.data_dict['SubjID'][good_subj_idx]
        self.data = data
        self.data_mat = data_mat
        self.good_keys = good_keys

        if not self.whiten:
            data_mat = (data_mat - data_mat.mean(0)) / data_mat.std(0)
        self.pca = PCA(whiten=self.whiten)
        self.pca.fit(np.abs(data_mat).T)

    def get_components(self):
        selected_idx = self.pca.explained_variance_ratio_ >= self.pc_threshhold
        return self.pca.components_[selected_idx]

    def get_projections(self):
        return np.dot(self.get_components(), self.data_mat)

    def report_asymmetry_loadings(self):
        mean_asymmetry = self.data_mat.std(axis=1)**2

        for pi, pc in enumerate(self.get_components()):
            print("%2d: (%.2f):" % (pi, self.pca.explained_variance_ratio_[pi]))
            sort_idx = np.argsort(np.abs(pc))
            for key, coeff, mag in zip(self.good_keys[sort_idx], pc[sort_idx], mean_asymmetry[sort_idx]):
                # if np.sign(mag) == hemi_sign:
                print("\t%s%.4f %-50s (%.3e / %.3e)" % (
                    ' ' if mag*coeff >= 0 else '',
                    coeff,
                    key,
                    coeff,
                    mag))

        return self.pca

    def report_behavior_correlations(self):
        # Now, take this factor and compare it to behavioral data
        # (or, should I throw the behavior in the PCA?)
        tbx_data = self.data.get_tbx_data()
        print('Found %d TBX keys.' % len(tbx_data))

        for pi, proj in enumerate(self.get_projections()):
            for ti, (key, val) in enumerate(tbx_data.items()):
                cur_idx = np.logical_not(np.isnan(val))
                rval, pval = pearsonr(proj[cur_idx], val[cur_idx])
                if pval < 0.05:
                    print('%d: %-20s (n=%4d) r=%.4f, p=%.4f' % (pi, key, cur_idx.sum(), rval, pval))

    def report_background_correlations(self):
        # Now, take this factor and compare it to hisory data
        # (or, should I throw the behavior in the PCA?)
        fdh_data = self.data.get_fdh_data()
        print('Found %d FDH keys.' % len(fdh_data))

        for pi, pc in enumerate(self.get_projections()):
            for ti, (key, val) in enumerate(fdh_data.items()):
                cur_idx = np.logical_not(np.isnan(val))
                rval, pval = pearsonr(proj[cur_idx], val[cur_idx])
                if pval < 0.05:
                    print('%d: %-20s (n=%4d) r=%.4f, p=%.4f' % (pi, key, cur_idx.sum(), rval, pval))
