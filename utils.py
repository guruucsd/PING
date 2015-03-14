import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


def asymmetry_index(left, right):
    """ Left and right should be arrays"""
    left = np.asarray(left)
    right = np.asarray(right)

    aidx = (left - right) / (left + right)
    aidx[np.isnan(aidx)] = 0
    return aidx


def do_and_plot_regression(X, Y, xlabel='', ylabel='', title='', ax=None):
    assert not np.any(np.isnan(X))
    assert not np.any(np.isnan(Y))

    m, b, rval,  pval, stderr = linregress(X, Y)
    if not ax:
        fh = plt.figure()
        ax = ax
    xlim = [X.size and X.min() or 0, X.size and X.max() or 1]
    xvals = np.linspace(xlim[0], xlim[1], 25)
    ax.plot(xvals, m * xvals + b, 'r', linewidth=3.)
    ax.set_title('%s\n(r=%.3f, p=%.3e, n=%d)' % (title, rval, pval, X.size))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.hold('on')
    ax.plot(xvals, xvals * 0., 'k--', linewidth=5)
    ax.scatter(X, Y)
    ax.set_ylim([-1, 1])
