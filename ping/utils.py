"""
Shared functions across scripts.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


def filter_dict(d, filter_fn):
    out_dict = dict()
    for key, val in d.items():
        if filter_fn(key, val):
            out_dict[key] = val

    return out_dict


def do_and_plot_regression(X, Y, covariates=[], xlabel='', ylabel='', title='', ax=None):
    assert not np.any(np.isnan(X))
    assert not np.any(np.isnan(Y))

    # Regress out stuff
    if covariates:
        raise NotImplementedException('covariates')
    m, b, rval,  pval, stderr = linregress(X, Y)
    if not ax:
        fh = plt.figure()
        ax = ax

    # Add standard deviation
    w_sz = 2.0
    xvals = np.arange(np.min(X) + w_sz / 2, np.max(X) - w_sz / 2, 0.01)
    yvals_mean = np.empty(xvals.shape)  # m * xvals + b
    yvals_std = np.empty(xvals.shape)
    for xi, xval in enumerate(xvals):
        idx = np.logical_and(xval - w_sz / 2 <= X, X <= xval + w_sz / 2)
        yvals_mean[xi] = Y[idx].mean()
        yvals_std[xi] = Y[idx].std()

    ax.plot(xvals, xvals * 0., 'k--', linewidth=5)
    ax.hold('on')

    if len(X) > 200:
        # ax.plot(xvals, yvals_mean, 'r', linewidth=3.)
        ax.fill_between(xvals, yvals_mean+yvals_std, yvals_mean-yvals_std,
                        facecolor='red', alpha=0.6)
    ax.scatter(X, Y)
    ax.plot(xvals, m * xvals + b, 'r', linewidth=3.)

    ax.set_title('%s\n(r=%.3f, p=%.3e, n=%d)' % (title, rval, pval, X.size))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim([-0.4, 0.4])

    return yvals_mean, yvals_std
