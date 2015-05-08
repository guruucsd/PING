"""
Functions related to asymmetry
"""
import numpy as np

from .access import (get_bilateral_hemi_keys, get_twohemi_keys,
                     get_nonhemi_prop_name)


def get_ai_prop_name(rh_prop_name):
    return '%s_AI' % get_nonhemi_prop_name(rh_prop_name)


def is_ai_prop_name(prop_name):
    return prop_name is not None and prop_name.endswith('_AI')


def asymmetry_index(left, right, mask_nan=True):
    """ Left and right should be arrays"""
    left = np.asarray(left)
    right = np.asarray(right)

    aidx = (left - right) / (left + right)
    if mask_nan:
        aidx[np.isnan(aidx)] = 0
    return aidx


def get_asymmetry_index(data, prop_name, mask_nan=True):
    """ Get the correponding left and right values for the prop_name,
    and returns the asymmetry index."""

    left_prop_name, right_prop_name = get_bilateral_hemi_keys(data, prop_name)

    # Select data within the group
    LH_data = np.asarray(data[left_prop_name].tolist())
    RH_data = np.asarray(data[right_prop_name].tolist())

    # Compute an asymmetry index
    prop_asymmetry = asymmetry_index(LH_data, RH_data, mask_nan=mask_nan)
    return prop_asymmetry


def compute_all_asymmetries(prefix):
    """ Loop over all properties to show asymmetry."""
    data = load_PING_data(scrub_fields=False)

    export_data = dict((('SubjID', data['SubjID'],),))

    # Process & plot the data.
    for prop_name in get_twohemi_keys(prefix, data.keys()):
        dest_prop_name = get_ai_prop_name(prop_name)
        export_data[dest_prop_name] = get_asymmetry_index(data, prop_name,
                                                          mask_nan=False)

    # Add one property for total asymmetry
    n_subj = len(data['SubjID'])
    for p in prefix:
        total_asymmetry = np.zeros((n_subj,))
        good_keys = filter(lambda k: k.startswith(p) and is_ai_prop_name(k),
                           export_data.keys())
        for key in good_keys:
            values = export_data[key].copy()
            values[np.isnan(values)] = 0.
            total_asymmetry += export_data[key]**2
        total_ai_key = get_ai_prop_name(p + '_TOTAL')
        export_data[total_ai_key] = np.sqrt(total_asymmetry)
        assert len(export_data[total_ai_key]) == n_subj

    return export_data
