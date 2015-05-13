"""
Plotting utilities
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_symmetric_matrix_as_triangle(mat, ax=None, labels=None, class_labels=None, vmin=0, vmax=1):
    """Plot symmetric matrix (like a covariance matrix) as a lower triangle."""

    # Scrub inputs
    if ax is None:
        ax = plt.figure().gca()
    if len(mat.shape) == 1:
        # Convert vector to matrix
        mat = scipy.spatial.distance.squareform(mat)
    mat[np.eye(mat.shape[0], dtype=bool)] = 0

    # Mask the matrix; will lead to transparency in imshow
    masked_mat = np.ma.masked_where(np.triu(np.eye(mat.shape[0])), mat)

    # interpolation: none needed for transparency...
    ax.set_axis_bgcolor(ax.get_figure().get_facecolor())
    img = ax.imshow(mat, vmin=vmin, vmax=vmax, interpolation='nearest')
    ax.set_frame_on(False)
    ax.tick_params(labelsize=16)

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.10)
    plt.colorbar(img, cax=cax)

    # Now label.
    if labels is None:
        ax.set_xticks([])
        ax.set_yticks([])

    elif class_labels is None or len(np.unique(class_labels)) == 1:
        sz = mat.shape[0]
        ax.set_xticks(range(sz))
        ax.set_xticklabels(class_labels, rotation='vertical')
        ax.set_yticks(range(sz))
        ax.set_yticklabels(class_labels)

    else:
        border_idx = [0]
        border_lbls = [class_labels[0]]
        for li, lbl_cls in enumerate(class_labels[1:]):
            if lbl_cls != border_lbls[-1]:
                print('Adding %s after %s' % (lbl_cls,
                                              border_lbls[-1]))
                border_idx.append(li)
                border_lbls.append(lbl_cls)
        ax.set_xticks(border_idx)
        ax.set_xticklabels(border_lbls, rotation='vertical')
        ax.set_yticks(border_idx)
        ax.set_yticklabels(border_lbls)



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
