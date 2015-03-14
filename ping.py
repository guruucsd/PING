"""
"""
import os

import numpy as np
import pandas

from utils import asymmetry_index


def col2prop(col_name):
    return col_name.replace('-', '.')


def load_PING_data():
    # Load data
    print("Loading data...")
    script_dir = os.path.abspath(os.path.dirname(__file__))
    csv_path = os.path.join(script_dir, 'PING_raw_data.csv')
    data = pandas.read_csv(csv_path)
    return data


def get_twohemi_keys(prefix, all_keys):
    """Given a key prefix, get all keys that have left/right pairs."""

    # Massage inputs
    if not isinstance(prefix, list):
        prefix = [prefix]
    prefix = np.asarray(prefix)
    rh_markers = ['Right', '.rh.', '.R_']

    good_keys = []
    for prop_name in all_keys:
        if (not np.any(np.asarray([r in prop_name for r in rh_markers])) or
                not np.any(np.asarray([p in prop_name for p in prefix]))):
            continue
        if 'vent' in prop_name.lower() and 'ventral' not in prop_name.lower():
            print("Skipping %s" % prop_name)
            continue
        good_keys.append(prop_name)
    return np.asarray(good_keys)



def get_twohemi_prop_names(data, prop_name):
    """ Given one hemisphere's property name,
    return both."""

    if 'Left' in prop_name or 'Right' in prop_name:
        left_prop_name = prop_name.replace('Right', 'Left')
        right_prop_name = prop_name.replace('Left', 'Right')
    elif '.lh.' in prop_name or '.rh.' in prop_name:
        left_prop_name = prop_name.replace('.rh.', '.lh.')
        right_prop_name = prop_name.replace('.lh.', '.rh.')
    elif '.L_' in prop_name or '.R_' in prop_name:
        left_prop_name = prop_name.replace('.R_', '.L_')
        right_prop_name = prop_name.replace('.L_', '.R_')
    else:
        raise ValueError("Unknown format for prop_name='%s'" % prop_name)

    return left_prop_name, right_prop_name


def get_asymmetry_index(data, prop_name):
    """ Get the correponding left and right values for the prop_name,
    and returns the asymmetry index."""

    left_prop_name, right_prop_name = get_twohemi_prop_names(data, prop_name)

    # Select data within the group
    LH_data = np.asarray(data[left_prop_name].tolist())
    RH_data = np.asarray(data[right_prop_name].tolist())

    # Compute an asymmetry index
    prop_asymmetry = asymmetry_index(LH_data, RH_data)
    return prop_asymmetry
