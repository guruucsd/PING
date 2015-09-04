"""
Various scatter plots

Goal is to have:
2D: typical scatter
3D: z is the size of the marker
4D: add in the color of the marker.

Should take ordered parameters for data keys on the input,
function should take keyword args.
"""
from argparse import ArgumentParser

import numpy as np
from matplotlib import pyplot as plt
from six import string_types

from ping.analysis.similarity import is_bad_key
from research.data import get_all_data, keytype2label


def parse_scatter_key(key):
    return key.split(':')


def compute_scatter_label(key, part=None):
    if key is None:
        return
    elif isinstance(key, string_types):
        key_type = keytype2label(parse_scatter_key(key)[0])
        method = parse_scatter_key(key)[1]
        if part is None:
           return '%s (%s)' % (key_type, method)
        elif part == 'key_type':
            return key_type
        elif part == 'method':
            return method
        else:
            raise NotImplementedError("Unrecognized part: %s" % part)  
    elif len(key) == 1:
        return compute_scatter_label(key[0], part=part)
    else:
        return [compute_scatter_label(k, part=part) for k in key]


def compute_key_data(data, key):
    # Now the hard part, ... interpreting the keys.
    # keys can be direct data arrays, *or* they can *across* keys (with operations).
    if key in data:
        # Easy case: it's just a data request!
        return data[key]
    elif len(parse_scatter_key(key)) != 2:
        # OK, it's nothing we know, error
        raise ValueError("key %s not found, nor of a computable format (suffix:function)" % key)
    else:
        # a key/operator pair
        suffix, op = parse_scatter_key(key)
        keys = [k for k in data if k.endswith(suffix)]
        assert len(keys) > 0, "Must find keys with filter!"
        f_data = {k: data[k][~np.isnan(data[k])] for k in keys}
        assert not np.any([np.any(np.isnan(v)) for v in f_data.values()]), "NO nan ANYWHERE..."
        return {k: getattr(f_data[k], op)() for k in keys}


def plot_scatter_4D(data, x_key, y_key, size_key=None, color_key=None, 
                    x_label=None, y_label=None, size_label=None, color_fn=None,
                    add_marker_text=False, ax=None):
    """y_key can be a list..."""
    colors = np.asarray(['b','r','g','y'])

    # Massage inputs
    if isinstance(y_key, string_types):
        y_key = [y_key]
    if isinstance(size_key, string_types):
        size_key = [size_key]
    if isinstance(color_key, string_types):
        color_key = [color_key]
    if x_label is None:
        x_label = compute_scatter_label(x_key)
    if y_label is None:
        y_label = compute_scatter_label(y_key)
    if size_label is None:
        size_label = compute_scatter_label(size_key)
    if ax is None:
        ax = plt.figure(figsize=(11, 10.5)).gca()

    # Now get all the data, and manipulate as needed
    kwargs = {
        'x': compute_key_data(data.data_dict, x_key),
        'y': [compute_key_data(data.data_dict, k) for k in y_key]}

    if size_key is not None:
        kwargs['s'] = np.asarray([compute_key_data(data.data_dict, k)
                                  for k in size_key])

    if color_key is not None:
        kwargs['c'] = np.asarray([compute_key_data(data.data_dict, k)
                                  for k in color_key])
    elif color_fn is not None:
        kwargs['c'] = np.asarray([{k: color_fn(k, v) for k, v in kwargs['x'].items()}])

    # Make sure everybody has the same keys
    common_keys = [data.get_nonhemi_key(k)
                   for k in kwargs['x'].keys()
                   if not is_bad_key(k) and ~np.all(np.isnan(kwargs['x'][k]))]
    if len(common_keys) == 0:
        raise ValueError('Your x key has an issue.')
    for key in list(set(kwargs.keys()) - set(['x'])):
        # Loop over
        cur_keys = [data.get_nonhemi_key(k)
                    for ddata in kwargs[key]
                    for k in ddata.keys()
                    if ~np.all(np.isnan(ddata[k]))]
        common_keys = [k for k in common_keys if k in cur_keys]
    if len(common_keys) == 0:
        raise ValueError('Your x and y keys have no overlap.')

    # Finally, we're safe to convert all of the data to numpy arrays,
    #   then massage the data.
    # NOTE: Loop over common keys, so all are ordered similarly
    #  BUT the actual keys in each dict is NOT the common_key,
    #  but some measure-specific version of it.
    #
    gmc = data.get_measure_key
    kwargs['x'] = np.asarray([kwargs['x'][gmc(ck, kwargs['x'].keys())]
                              for ck in common_keys])
    for key in list(set(kwargs.keys()) - set(['x'])):
        kwargs[key] = np.asarray([sdata[gmc(ck, sdata.keys())]
                                  for sdata in kwargs[key]
                                  for ck in common_keys])

    if 's' in kwargs:
        kwargs['s'] = 1000 * kwargs['s'] / np.abs(kwargs['s']).mean()
    if 'c' in kwargs:
        kwargs['c'] = colors[kwargs['c']].ravel()

    # Now plot it, and annotate it!
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

        for label, x, y, s in zip(common_keys, kwargs['x'], kwargs['y'], kwargs['s']):
            locs = locals()
            annotations = [key for key, sval in zip(['x', 'y', 's'], [1.35, 1.5, 2])
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
    return ax


def do_scatter(prefix, x_key, y_key, size_key=None, color_key=None, dataset='ping'):

    prefix = prefix.split(',')
    y_key = y_key.split(',')
    size_key = size_key.split(',') if size_key else size_key
    color_key = color_key.split(',') if color_key else color_key

    # Get prefix
    prefix_filter_fn = lambda k, v: np.any([k.startswith(p) for p in prefix])

    # Load the data (should group, but ... later.)
    data = get_all_data(dataset).filter(prefix_filter_fn)

    if size_key is not None:
        size_label = ' Marker size indicates\n %s %s' % (
            compute_scatter_label(size_key, part='key_type').lower(),
             ', '.join([data.prefix2text(p).lower() for p in prefix]))
    else:
        size_label = None
    ax = plot_scatter_4D(data, x_key=x_key, y_key=y_key, size_key=size_key, color_key=color_key,
                         size_label=size_label, add_marker_text=True)
    # x_label='Asymmetry Index (mean)', y_label='Asymmetry Index (std)',

    ax.get_figure().suptitle(', '.join([data.prefix2text(p)
                                        for p in prefix]),
                             fontsize=24)
    plt.show()


if __name__ == '__main__':
    axis_choices = ['AI:mean', 'AI:std',
                    'LH_PLUS_RH:mean', 'LH_PLUS_RH:std',
                    'TOTAL:mean', 'TOTAL:std']
    parser = ArgumentParser(description="Scatter plot on any two data"
                            " arrays, with additional data arrays that"
                            " optionally control marker size and color.")
    parser.add_argument('prefix', help="comma-separated list of prefixes to"
                                       " include in the analysis")
    parser.add_argument('x_key', choices=axis_choices)
    parser.add_argument('y_key', choices=axis_choices)
    parser.add_argument('size_key', choices=axis_choices,
                        nargs='?', default=None)
    parser.add_argument('color_key', choices=axis_choices,
                        nargs='?', default=None)
    parser.add_argument('--dataset', choices=['ping', 'destrieux'],
                        nargs='?', default='ping')
    args = parser.parse_args()
    do_scatter(**vars(args))
