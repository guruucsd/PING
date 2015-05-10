"""
Plotting utilities
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial


def plot_symmetric_matrix_as_triangle(mat, ax=None, lbls=None, vmin=0, vmax=1):
    """Plot symmetric matrix (like a covariance matrix) as a lower triangle."""

    # Scrub inputs
    if ax is None:
        ax = plt.figure().gca()
    if len(mat.shape) == 1:
        # Convert vector to matrix
        mat = scipy.spatial.distance.squareform(mat)

    # Mask the matrix; will lead to transparency in imshow
    masked_mat = np.ma.masked_where(np.triu(np.eye(mat.shape[0])), mat)

    # interpolation: none needed for transparency...
    ax.set_axis_bgcolor(ax.get_figure().get_facecolor())
    ax.imshow(mat, vmin=vmin, vmax=vmax, interpolation='none')
    ax.set_frame_on(False)

    if lbls:
        sz = mat.shape[0]
        ax.set_xticks(range(sz - 1))
        ax.set_xticklabels(lbls[:sz - 1])
        ax.set_yticks(range(1, sz))
        ax.set_yticklabels(lbls[1:])


def equalize_xlims(fh):
    xlims = np.asarray([ax.get_xlim() for ax in fh.get_axes()])
    xlim = [xlims[:, 0].min(), xlims[:, 1].max()]

    for ax in fh.get_axes():
        ax.set_xlim(xlim)


def equalize_ylims(fh):
    ylims = np.asarray([ax.get_ylim() for ax in fh.get_axes()])
    ylim = [ylims[:, 0].min(), ylims[:, 1].max()]

    for ax in fh.get_axes():
        ax.set_ylim(ylim)
