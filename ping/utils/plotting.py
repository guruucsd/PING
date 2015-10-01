"""
Plotting utilities
"""
import matplotlib.patches as pat
import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot_normalized_hist(data, ax=None, **kwargs):
    ax = ax or plt.figure().gca()
    
    x, b, h = ax.hist(data, normed=True, **kwargs)
    if len(h) == 0:
        pass
    elif isinstance(h[0], pat.Rectangle):
        for item in h:
            item.set_height(item.get_height() / sum(x))
    else:
        for item, dat in zip(h, x):
            for rect in item.patches:
                rect.set_height(rect.get_height() / sum(dat))

    ax.set_ylim([0, 1])

    return x, b, h


def plot_symmetric_matrix_as_triangle(mat, ax=None, labels=None,
                                      class_labels=None, vmin=0, vmax=1,
                                      plotengine='matplotlib'):
    """Plot symmetric matrix (like a covariance matrix) as a lower triangle.

    Can accept matrix in vector or matrix form."""

    ## Prep data and labels

    if len(mat.shape) == 1:
        # Convert vector to matrix
        mat = scipy.spatial.distance.squareform(mat)
    # mat[np.eye(mat.shape[0], dtype=bool)] = 0

    # Mask the matrix; will lead to transparency in imshow
    masked_mat = np.ma.masked_where(np.tril(np.ones(mat.shape)), mat)
    masked_mat = masked_mat.T #[::-1, ::-1]

    ## Now plot
    ticks = dict()

    # Now label.
    if labels is None:
        ticks = dict(x=[], y=[], xlab=[], ylab=[])

    elif class_labels is None or len(np.unique(class_labels)) == 1:
        sz = mat.shape[0]
        ticks = dict(x=range(sz - 1), xlab=labels[:-1],
                     y=range(1, sz), ylab=labels[1:])

    else:
        border_idx = [0]
        border_lbls = [class_labels[0]]
        for li, lbl_cls in enumerate(class_labels[1:]):
            if lbl_cls != border_lbls[-1]:
                print('Adding %s after %s' % (lbl_cls,
                                              border_lbls[-1]))
                border_idx.append(li)
                border_lbls.append(lbl_cls)
        ticks = dict(x=border_idx, xlab=border_lbls,
                     y=border_idx, ylab=border_lbls)



    ## Plotting
    if plotengine in ['matplotlib', 'mpld3']:
        # Scrub inputs
        if ax is None:
            ax = plt.figure().gca()

        # interpolation: none needed for transparency...
        ax.set_axis_bgcolor(ax.get_figure().get_facecolor())
        img = ax.imshow(masked_mat, vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.set_frame_on(False)
        ax.tick_params(labelsize=16)

        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.10)
        plt.colorbar(img, cax=cax)

        ax.set_xticks(ticks['x'])
        ax.set_yticks(ticks['y'])
        ax.set_xticklabels(ticks['xlab'])
        ax.set_yticklabels(ticks['ylab'])
        ax.xaxis.set_ticks_position('none')
        ax.yaxis.set_ticks_position('none')

    elif plotengine in ['bokeh']:
        from collections import OrderedDict
        from bokeh.plotting import figure
        from bokeh.models import HoverTool, ColumnDataSource
        from bokeh.sampledata.les_mis import data

        #nodes = data['nodes']
        # names = [node['name'] for node in sorted(data['nodes'], key=lambda x: x['group'])]

        # N = len(nodes)
        # counts = np.zeros((N, N))
        # for link in data['links']:
        #     counts[link['source'], link['target']] = link['value']
        #     counts[link['target'], link['source']] = link['value']

        colormap = [
            "#444444", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
            "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"
        ]
        # import pdb; pdb.set_trace()
        xname = []
        yname = []
        color = []
        alpha = []
        for i, n1 in enumerate(labels):
            for j, n2 in enumerate(labels):
                xname.append(n1)
                yname.append(n2)

                a_offset = 0.025
                val = (1. - a_offset)* abs(mat[i,j]) / max(abs(vmax), abs(vmin))  # normalize to 0..0.9 range
                a = min(val, (1. - a_offset)) + a_offset
                alpha.append(a)

                if class_labels[i] == class_labels[j]:
                    unique_class_labels = np.unique(class_labels).tolist()
                    class_idx = unique_class_labels.index(class_labels[i])
                    # color.append(colormap[1 + class_idx])
                    color.append('#f00' if mat[i,j] >= 0 else '#00f')
                else:
                    color.append(colormap[0])

        source = ColumnDataSource(
            data=dict(
                xname=xname,
                yname=yname,
                colors=color,
                alphas=alpha,
                count=mat.flatten()))

        p = figure(x_axis_location="above", tools="resize,hover,save",
                   x_range=list(reversed(labels)), y_range=labels)
        p.plot_width = 800
        p.plot_height = 800

        p.rect('xname', 'yname', 0.9, 0.9, source=source,
               color='colors', alpha='alphas', line_color=None)

        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.major_label_text_font_size = "5pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = np.pi/3

        hover = p.select(dict(type=HoverTool))
        hover.tooltips = OrderedDict([('names', '@yname, @xname')])

        ax = p
    return ax


def equalize_xlims(fh, xlim=None):
    if xlim is None:
        xlims = np.asarray([ax.get_xlim() for ax in fh.get_axes()])
        xlim = [xlims[:, 0].min(), xlims[:, 1].max()]

    for ax in fh.get_axes():
        ax.set_xlim(xlim)


def equalize_ylims(fh, ylim=None):
    if ylim is None:
        ylims = np.asarray([ax.get_ylim() for ax in fh.get_axes()])
        ylim = [ylims[:, 0].min(), ylims[:, 1].max()]

    for ax in fh.get_axes():
        ax.set_ylim(ylim)
