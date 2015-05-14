"""
Various scatter plots

Goal is to have:
2D: typical scatter
3D: z is the size of the marker
4D: add in the color of the marker.

Should take ordered parameters for data keys on the input,
function should take keyword args.
"""
import copy
import sys
from collections import OrderedDict
from functools import partial

import numpy as np
from matplotlib import pyplot as plt
from six import string_types

from ping.analysis.similarity import (is_bad_key,
                                      compare_similarity_matrices,
                                      compute_similarity_matrices,
                                      visualize_similarity_matrices)
from ping.data import (PINGData, which_hemi, get_nonhemi_key, is_nonimaging_key, 
                       get_anatomical_name, anatomical_sort, get_prefixes)
from ping.utils import filter_dict
from research.asymmetry import is_ai_key
from research.data import get_all_data
from research.grouping import parse_filter_args


def compute_key_data(data, key):
    # Now the hard part, ... interpreting the keys.
    # keys can be direct data arrays, *or* they can *across* keys (with operations).
    if key in data:
        # Easy case: it's just a data request!
        return data[key]
    elif ':' not in key:
        # OK, it's nothing we know, error
        raise ValueError("key %s not found, nor of a computable format (suffix:function)" % key)
    else:
        # a key/operator pair
        suffix, op = key.split(':')
        keys = [k for k in data if k.endswith(suffix)]
        f_data = {k: data[k][~np.isnan(data[k])] for k in keys}
        assert not np.any([np.any(np.isnan(v)) for v in f_data.values()]), "NO nan ANYWHERE..."

        return {k: getattr(f_data[k], op)() for k in keys}


def plot_scatter_4D(data, x_key, y_key, size_key=None, color_key=None, 
                    x_labels=None, y_labels=None, color_fn=None,
                    ax=None):
    """y_key can be a list..."""
    colors = np.asarray(['b','r','g','y'])

    # Massage inputs
    if isinstance(y_key, string_types):
        y_key = [y_key]
    if isinstance(size_key, string_types):
        size_key = [size_key]
    if isinstance(color_key, string_types):
        color_key = [color_key]
    if ax is None:
        ax = plt.figure().gca()

    # Now get all the data, and manipulate as needed
    kwargs = {
        'x': compute_key_data(data.data_dict, x_key),
        'y': [compute_key_data(data.data_dict, k) for k in y_key]}

    if size_key:
        kwargs['s'] = np.asarray([compute_key_data(data.data_dict, k)
                                  for k in size_key])

    if color_key:
        kwargs['c'] = np.asarray([compute_key_data(data.data_dict, k) 
                                  for k in color_key])
    elif color_fn:
        kwargs['c'] = np.asarray([{k: color_fn(k, v) for k, v in kwargs['x'].items()}])

    # Make sure everybody has the same keys
    common_keys = [get_nonhemi_key(k) 
                   for k in kwargs['x'].keys()
                   if not is_bad_key(k) and ~np.isnan(kwargs['x'][k])]

    for key in list(set(kwargs.keys()) - set(['x'])):
        # Loop over 
        cur_keys = [get_nonhemi_key(k)
                    for ddata in kwargs[key]
                    for k in ddata.keys()
                    if ~np.isnan(ddata[k])]
        common_keys = [k for k in common_keys if k in cur_keys]
    print(common_keys)

    # Finally, we're safe to convert all of the data to numpy arrays,
    #   then massage the data.
    kwargs['x'] = np.asarray([val for k, val in kwargs['x'].items()
                              if np.any([ck in k for ck in common_keys])])
    for key in list(set(kwargs.keys()) - set(['x'])):
        kwargs[key] = np.asarray([val for sdata in kwargs[key]
                                      for k, val in sdata.items()
                                      if np.any([ck in k for ck in common_keys])])

    kwargs['s'] = 1000  * kwargs['s'] / np.abs(kwargs['s']).mean()
    kwargs['c'] = colors[kwargs['c']].ravel()

    # Now plot it, and annotate it!
    ax.scatter(**kwargs)

    # Interesting if it's outside of some range of values
    is_interesting = lambda v, varr, dist=1.5: np.abs(varr.mean() - v) >= dist * varr.std()

    for label, x, y, s, c in zip(common_keys, kwargs['x'], kwargs['y'], kwargs['s'], kwargs['c']):
        annotations = [key for key, sval in zip(['x', 'y', 's'], [1.5, 1.5, 2])
                       if is_interesting(locals()[key], kwargs[key], sval)]
        if len(annotations) > 0:
            plt.annotate(
                '%s (%s)' % (get_anatomical_name(get_nonhemi_key(label)), ', '.join(annotations)),
                xy = (x, y), xytext = (-20, 20),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=0.5', fc = c, alpha = 0.5),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'),
                )

    plt.show()

    return ax


def do_usage(args, error_msg=None):
    if error_msg is not None:
        print("*** ERROR *** : %s" % error_msg)
    print("\nUsage: %s prefix x_key y_key [size_key] [color_key]" % args[0])
    print("\tScatter plot on any two data arrays, with additional data arrays that")
    print("\toptionally control marker size and color.")
    print("\n\tprefix: asdfadf.")
    print("\tx_key: data key to control x values.")
    print("\ty_key: comma-separated list of keys for y series.")
    print("\tsize_key: (optional) comma-separated list of keys for controlling size.")
    print("\t\tNote: must be one key, or as many keys as in y_key.")
    print("\tcolor_key: (optional) comma-separated list of keys for controlling color.")
    print("\t\tNote: must be one key, or as many keys as in y_key.")


if __name__ != '__main__':
    pass

elif len(sys.argv) > 6:
    do_usage("Too many arguments.")

elif len(sys.argv) < 4:
    do_usage("Too few keys.")

else:

    prefix, x_key, y_key = sys.argv[1:4]
    prefix = prefix.split(',')
    y_key = y_key.split(',')
    size_key = None if len(sys.argv) < 5 else sys.argv[4].split(',')
    color_key = None if len(sys.argv) < 6 else sys.argv[5].split(',')

    # Get prefix
    prefix_filter_fn = lambda k, v: np.any([k.startswith(p) for p in prefix])

    # Load the data (should group, but ... later.)
    data = get_all_data().filter(prefix_filter_fn)

    plot_scatter_4D(data, x_key=x_key, y_key=y_key, size_key=size_key, color_key=color_key,
        color_fn=lambda k, v: np.nonzero([k.startswith(p) for p in prefix])[0]) 
