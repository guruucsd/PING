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


def do_and_plot_regression(X, Y, covariates=[], xlabel='', ylabel='',
                           title='', ax=None):
    assert not np.any(np.isnan(X))
    assert not np.any(np.isnan(Y))

    # Regress out stuff
    if covariates:
        raise NotImplementedException('covariates')
    m, b, rval, pval, stderr = linregress(X, Y)
    if not ax:
        fh = plt.figure()
        ax = ax

    # Add standard deviation
    w_sz = 3.0
    xvals = np.arange(np.min(X) + w_sz / 2, np.max(X) - w_sz / 2, 0.01)
    yvals_mean = np.empty(xvals.shape)  # m * xvals + b
    yvals_std = np.empty(xvals.shape)
    for xi, xval in enumerate(xvals):
        idx = np.logical_and(xval - w_sz / 2 <= X, X <= xval + w_sz / 2)
        yvals_mean[xi] = Y[idx].mean()
        yvals_std[xi] = Y[idx].std()

    ax.plot([2, 22], [0., 0], 'k--', linewidth=5)
    ax.hold('on')

    if len(X) > 200:
        # ax.plot(xvals, yvals_mean, 'r', linewidth=3.)
        ax.fill_between(xvals, yvals_mean+yvals_std, yvals_mean-yvals_std,
                        facecolor=[0., 0., 0., 0.4])
    ax.scatter(X, Y)
    linvals = np.asarray([X.min(), X.max()])
    ax.plot(linvals, m * linvals + b, c=[1, 0., 0., 0.8], linewidth=7.)

    ax.set_title('%s\n(r=%.3f, p=%.3f)' % (title, rval, pval),
                 fontsize=24)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)

    ax.tick_params(labelsize=16)
    ax.set_ylim([-0.4, 0.4])
    ax.set_xlim([2, 22])

    return m, b, rval, pval
