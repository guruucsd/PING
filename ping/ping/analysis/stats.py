"""
"""
import numpy as np
import scipy.stats


def pdist_nan(data_mat, metric='correlation', standardize=False):
    if not np.isnan(data_mat.sum()):
        return scipy.spatial.distance.pdist(data_mat, 'correlation')

    else:
        vec = []
        n_keys = data_mat.shape[0]
        for ki in range(n_keys):
            for li in range(ki + 1, n_keys):
                good_idx = ~np.isnan(data_mat[ki] + data_mat[li])
                kvals = data_mat[ki][good_idx]
                lvals = data_mat[li][good_idx]
                if standardize:
                    kvals = scipy.stats.mstats.zscore(kvals)
                    lvals = scipy.stats.mstats.zscore(lvals)

                vec.append(scipy.stats.pearsonr(kvals, lvals)[0])
                if np.isnan(vec[-1]):
                    # assert kvals.std() == 0 or lvals.std() == 0
                    vec[-1] = 0
                assert np.abs(vec[-1]) <= 1
        return 1 - np.asarray(vec)
