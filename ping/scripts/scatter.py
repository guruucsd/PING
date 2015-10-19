"""
Various scatter plots

Goal is to have:
2D: typical scatter
3D: z is the size of the marker
4D: add in the color of the marker.

Should take ordered parameters for data keys on the input,
function should take keyword args.
"""
import simplejson
from collections import OrderedDict

import numpy as np
from matplotlib import pyplot as plt
from six import string_types

from ..ping.analysis.similarity import is_bad_key
from ..ping.apps import PINGSession
from ..ping.data import prefix2measure
from ..research.apps import ResearchArgParser
from ..research.data import get_all_data, keytype2label
from ..research.plotting import show_plots


def parse_scatter_key(key):
    parts = key.split(':')
    if len(parts) == 2:
        return dict(zip(('suffix', 'operator'), parts))
    elif len(parts) == 3:
        d = dict(zip(('prefix', 'suffix', 'operator'), parts))
        d['quantity'] = prefix2measure(d['prefix'])
        return d
    else:
        raise Exception("Unexpected key format: %s" % key)


def compute_scatter_label(key, part=None):
    # part: part of key to return
    if key is None:
        return None

    if isinstance(key, string_types):
        parts = parse_scatter_key(key)
        quantity = parts.get('quantity', None)
        key_type = keytype2label(parts['suffix'])
        method = parse_scatter_key(key)['operator']

        if part == 'quantity':
            return quantity
        if part == 'key_type':
            return key_type
        if part == 'method':
            return method
        if part is None and quantity:
           return '%s (%s %s)' % (quantity, method, key_type)
        if part is None:
           return '%s %s' % (method, key_type)
        raise NotImplementedError("Unrecognized part: %s" % part)

    # Lists
    if len(key) == 1:
        return compute_scatter_label(key[0], part=part)
    return [compute_scatter_label(k, part=part) for k in key]


def compute_key_data(data, key):
    # Now the hard part, ... interpreting the keys.
    # keys can be direct data arrays, *or* they can *across* keys (with operations).
    if key in data:
        # Easy case: it's just a data request!
        return data[key]

    parts = parse_scatter_key(key)
    new_keys = dict([(key, key) for key in data])

    if 'prefix' in parts:
        new_keys = dict([(k, nk[len(parts['prefix']):])
                         for k, nk in new_keys.items()
                        if k.startswith(parts['prefix'])])
        assert len(new_keys) > 0, "Must find keys with filter!"

    if 'suffix' in parts:
        new_keys = dict([(k, nk[:-len(parts['suffix'])])
                         for k, nk in new_keys.items()
                         if k.endswith(parts['suffix'])])
        assert len(new_keys) > 0, "Must find keys with filter!"

    # eliminate keys with nan data
    f_data = {k: data[k][~np.isnan(data[k])] for k in new_keys}  # remove nan subjects
    assert not np.any([np.any(np.isnan(v)) for v in f_data.values()]), "NO nan ANYWHERE..."

    return OrderedDict([(nk, getattr(f_data[k], parts['operator'])())
                        for k, nk in new_keys.items()])


def decimate_data(data, x_key, y_key, size_key=None, color_key=None, color_fn=None):

    # Massage inputs
    if isinstance(y_key, string_types):
        y_key = [y_key]
    if isinstance(size_key, string_types):
        size_key = [size_key]
    if isinstance(color_key, string_types):
        color_key = [color_key]

    # Now get all the data, and manipulate as needed
    out_data = {
        'x': compute_key_data(data.data_dict, x_key),
        'y': [compute_key_data(data.data_dict, k) for k in y_key]}

    if size_key is not None:
        out_data['s'] = np.asarray([compute_key_data(data.data_dict, k)
                                  for k in size_key])

    if color_key is not None:
        out_data['c'] = np.asarray([compute_key_data(data.data_dict, k)
                                  for k in color_key])
    elif color_fn is not None:
        out_data['c'] = np.asarray([{k: color_fn(k, v) for k, v in out_data['x'].items()}])

    # Make sure everybody has the same keys
    common_keys = [data.get_nonhemi_key(k)
                   for k in out_data['x'].keys()
                   if not is_bad_key(k) and ~np.all(np.isnan(out_data['x'][k]))]
    if len(common_keys) == 0:
        raise ValueError('Your x key has an issue.')
    for key in list(set(out_data.keys()) - set(['x'])):
        # Loop over
        cur_keys = [data.get_nonhemi_key(k)
                    for ddata in out_data[key]
                    for k in ddata.keys()
                    if ~np.all(np.isnan(ddata[k]))]
        common_keys = [k for k in common_keys if k in cur_keys]
    if len(common_keys) == 0:
        raise ValueError('Your keys (%s) have no overlap.' % out_data.keys())

    # Finally, we're safe to convert all of the data to numpy arrays,
    #   then massage the data.
    # NOTE: Loop over common keys, so all are ordered similarly
    #  BUT the actual keys in each dict is NOT the common_key,
    #  but some measure-specific version of it.
    #
    gmc = data.get_measure_key
    out_data['x'] = np.asarray([out_data['x'][gmc(ck, out_data['x'].keys())]
                              for ck in common_keys])
    for key in list(set(out_data.keys()) - set(['x'])):
        out_data[key] = np.asarray([sdata[gmc(ck, sdata.keys())]
                                  for sdata in out_data[key]
                                  for ck in common_keys])

    if 's' in out_data:
        out_data['s'] = 1000 * out_data['s'] / np.abs(out_data['s']).mean()
    if 'c' in out_data:
        out_data['c'] = colors[out_data['c']].ravel()

    out_data.update(dict(keys=common_keys))
    return out_data


def plot_scatter_4D(data, x_key, y_key, size_key=None, color_key=None,
                    x_label=None, y_label=None, size_label=None, color_fn=None,
                    add_marker_text=False, ax=None, title=None,
                    plotengine='matplotlib'):
    """y_key can be a list..."""
    colors = np.asarray(['b','r','g','y'])

    if x_label is None:
        x_label = compute_scatter_label(x_key)
    if y_label is None:
        y_label = compute_scatter_label(y_key)
    if size_label is None:
        size_label = compute_scatter_label(size_key)

    # Data is in format needed to scatter, so we call it
    #   kwargs.
    kwargs = decimate_data(data, x_key=x_key, y_key=y_key,
                           size_key=size_key, color_key=color_key,
                           color_fn=color_fn)
    common_keys = kwargs.pop('keys')

    # Now plot it, and annotate it!
    if plotengine in ['matplotlib', 'mpld3']:
        if ax is None:
            ax = plt.figure(figsize=(11, 10.5)).gca()
        ax.scatter(**kwargs)
        ax.tick_params(labelsize=16)
        if x_label:
            ax.set_xlabel(x_label, fontsize=18)
        if y_label:
            ax.set_ylabel(y_label, fontsize=18)
        if size_label:
            if 'thickness' in size_label:  # hack
                loc = 'upper left'
            else:
                loc = 'upper right'
            ax.legend([size_label], loc=loc)

        if add_marker_text:
            # Interesting if it's outside of some range of values
            is_interesting = lambda v, varr, dist: np.abs(varr.mean() - v) >= dist * varr.std()

            for xi, (label, x, y) in enumerate(zip(common_keys, kwargs['x'], kwargs['y'])): # , kwargs['s']):
                locs = locals()
                annotations = [key for key, sval in zip(['x', 'y'], [1.35, 1.5])
                               if is_interesting(locs[key], kwargs[key], sval)]
                if len(annotations) > 0:
                    plt.annotate(
                        '%s (%s)' % (data.get_anatomical_name(data.get_nonhemi_key(label)), ', '.join(annotations)),
                        xy=(x, y), xytext=(25, 25),
                        textcoords='offset points', ha='right', va='bottom',
                        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                        arrowprops=dict(arrowstyle='->', connectionstyle='arc3, rad=0'),
                        fontsize=16)

        plt.axis('equal') #ax.set_aspect('equal')
        if np.any(kwargs['x'] <= 0) and np.any(kwargs['x'] >= 0):
            ax.plot([0, 0], ax.get_ylim(), 'k--')  # y-axis
        if np.any(kwargs['y'] <= 0) and np.any(kwargs['y'] >= 0):
           ax.plot(ax.get_xlim(), [0, 0], 'k--')  # x-axis
        plt.axis('tight')

        ax.get_figure().suptitle(title, fontsize=24)
    elif plotengine in ['bokeh', 'bokeh-silent']:
        from bokeh.io import show
        from bokeh.models.glyphs import Circle
        from bokeh.models import (Plot, DataRange1d, LinearAxis, Legend,
                                  ColumnDataSource, PanTool, WheelZoomTool,
                                  HoverTool, CrosshairTool, PreviewSaveTool,
                                  ResizeTool, BoxZoomTool)
        source = ColumnDataSource(data=dict(
            label=[data.get_anatomical_name(data.get_nonhemi_key(label))
                   for label in common_keys],
            x=kwargs['x'],
            y=kwargs['y'],
            s=kwargs.get('s', np.ones(kwargs['x'].shape) * 1000) / 1.E5))
        xdr = DataRange1d()
        ydr = DataRange1d()
        plot = Plot(x_range=xdr, y_range=ydr, name='scatter', title=title)
        circle = Circle(x="x", y="y", radius="s",
                        fill_color=kwargs.get('c', 'blue'),
                        line_color="black")
        circle_renderer = plot.add_glyph(source, circle)

        if size_label:
            legend = Legend(orientation="top_right")
            legend.legends = [(size_label, [circle_renderer])]
            plot.add_layout(legend)

        plot.add_layout(LinearAxis(axis_label=x_label), 'below')
        plot.add_layout(LinearAxis(axis_label=y_label), 'left')
        plot.add_tools(PanTool(), WheelZoomTool(), CrosshairTool(),
                       PreviewSaveTool(), ResizeTool(), BoxZoomTool(),
                       HoverTool(tooltips=[('label', '@label')]))
        ax = plot
    return ax


def do_scatter(prefixes, x_key, y_key, size_key=None, color_key=None,
               atlas='desikan', username=None, passwd=None,
               output_format='matplotlib', data_dir='data', output_dir='data'):

    y_key = y_key.split(',')
    size_key = size_key.split(',') if size_key else size_key
    color_key = color_key.split(',') if color_key else color_key

    # Load the data (should group, but ... later.),
    # then filter by prefix
    data = get_all_data(atlas, username=username, passwd=passwd, data_dir=data_dir)
    data = data.filter(lambda k, v: np.any([k.startswith(p) for p in prefixes]))

    if size_key is not None:
        size_label = ' Marker size indicates\n %s %s' % (
            compute_scatter_label(size_key, part='key_type').lower(),
             ', '.join([data.prefix2text(p).lower() for p in prefixes]))
    else:
        size_label = None

    if output_format in ['json']:
        scatter_data = decimate_data(data, x_key=x_key, y_key=y_key, size_key=size_key,
                                     color_key=color_key)
        keys = [data.get_anatomical_name(data.get_nonhemi_key(key))
                for key in scatter_data.pop('keys')]
        out_dict = dict()
        for k, v in scatter_data.items():
            out_dict[k] = dict(zip(keys, v))
        out_file = '%s_scatter.json' % ','.join(prefixes)
        with open(out_file, 'wb') as fp:
            simplejson.dump(out_dict, fp)

    else:
        ax = plot_scatter_4D(data, x_key=x_key, y_key=y_key, size_key=size_key,
                             color_key=color_key, size_label=size_label,
                             add_marker_text=True,
                             title=', '.join([data.prefix2text(p)
                                              for p in prefixes]),
                             plotengine=output_format)
        # x_label='Asymmetry Index (mean)', y_label='Asymmetry Index (std)',

        show_plots(output_format, ax=ax, output_dir=output_dir)


if __name__ == '__main__':
    axis_choices = ['AI:mean', 'AI:std',
                    'LH_PLUS_RH:mean', 'LH_PLUS_RH:std']
    parser = ResearchArgParser(description="Scatter plot on any two data"
                               " arrays, with additional data arrays that"
                               " optionally control marker size and color.",
                               common_args=['prefixes',
                                            'atlas', 'username', 'passwd',
                                            'output-dir'])
    parser.add_argument('x_key')#, choices=ResearchArgParser.axis_choices)
    parser.add_argument('y_key')#, choices=ResearchArgParser.axis_choices)
    parser.add_argument('size_key',#, choices=ResearchArgParser.axis_choices,
                        nargs='?', default=None)
    parser.add_argument('color_key',# choices=ResearchArgParser.axis_choices,
                        nargs='?', default=None)
    parser.add_argument('--output-format', nargs='?', default='matplotlib',
                        choices=['matplotlib', 'mpld3', 'json',
                                 'bokeh', 'bokeh-silent'])
    args = parser.parse_args()
    do_scatter(**vars(args))
