"""
Shared functions across scripts.
"""
import copy

import numpy as np
import matplotlib.pyplot as plt
import seaborn
from scipy.stats import linregress


def filter_dict(d, filter_fn):
    out_dict = dict()
    for key, val in d.items():
        if filter_fn(key, val):
            out_dict[key] = val

    return out_dict


def do_and_plot_regression(X, Y, covariates=[], xlabel=None, ylabel=None,
                           title=None, ax=None, xlim=None, ylim=None,
                           colori=0, show_std=False, plotengine='matplotlib'):
    """
    Parameters:
        colori : (optional, default=0) which color to choose
            (0=blue scheme, 1=red scheme, 2=green scheme)
        show_std: (optional, default=False) show filled standard deviation?
    """
    assert not np.any(np.isnan(X))
    assert not np.any(np.isnan(Y))

    colors = np.asarray([[1., 0., 0., 1.], [0., 0., 1., 1.], [0., 1., 0., 1.]])

    # Regress out stuff
    if covariates:
        raise NotImplementedException('covariates')
    m, b, rval, pval, stderr = linregress(X, Y)
    if not ax:
        fh = plt.figure()
        ax = ax

    # Add standard deviation
    w_sz = (np.max(X) - np.min(X)) / 25.  # window size

    xvals = np.linspace(np.min(X) + w_sz / 2, np.max(X) - w_sz / 2, 1000)
    yvals_mean = np.empty(xvals.shape)  # m * xvals + b
    yvals_std = np.empty(xvals.shape)
    for xi, xval in enumerate(xvals):
        idx = np.logical_and(xval - w_sz / 2 <= X, X <= xval + w_sz / 2)
        if idx.sum() > 0:
            yvals_mean[xi] = Y[idx].mean()
            yvals_std[xi] = Y[idx].std()

    if plotengine in ['matplotlib', 'mpld3']:
        ax.plot([2, 22], [0., 0], 'k--', linewidth=5)  # axes
        ax.hold('on')

        if show_std:
            std_colors = copy.copy(colors)
            std_colors[:, -1] = 0.4  # alpha
            #std_colors[std_colors == 0.] = 0.4  # faded color
            ax.fill_between(xvals, yvals_mean+yvals_std, yvals_mean-yvals_std,
                            facecolor=std_colors[colori])

        ax.scatter(X, Y, c=colors[colori], s=35)

        reg_colors = copy.copy(colors)
        reg_colors[:, -1] = 0.8  # alpha
        linvals = np.asarray([X.min(), X.max()])
        ax.plot(linvals, m * linvals + b, c=reg_colors[colori], linewidth=7.)

        # add metadata
        if title:
            ax.set_title('%s\n(r=%.3f, p=%.3f)' % (title, rval, pval),
                         fontsize=24)
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=18)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=18)

        ax.tick_params(labelsize=16)
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
    else:
        raise NotImplementedException()

    return m, b, rval, pval
