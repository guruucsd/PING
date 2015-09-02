"""
Functions related to asymmetry
"""
import numpy as np

from ping.data import PINGData


def get_ai_key(rh_key):
    return '%s_AI' % PINGData.get_nonhemi_key(rh_key)


def is_ai_key(key):
    return key is not None and key.endswith('_AI')


def asymmetry_index(left, right, mask_nan=False):
    """ Left and right should be arrays"""
    left = np.asarray(left)
    right = np.asarray(right)

    aidx = (left - right) / (left + right)
    if mask_nan:
        aidx[np.isnan(aidx)] = 0
    return aidx


def get_asymmetry_index(data, key, mask_nan=False):
    """ Get the correponding left and right values for the key,
    and returns the asymmetry index."""

    left_key, right_key = PINGData.get_bilateral_hemi_keys(key)

    # Select data within the group
    LH_data = np.asarray(data[left_key].tolist())
    RH_data = np.asarray(data[right_key].tolist())

    # Compute an asymmetry index
    prop_asymmetry = asymmetry_index(LH_data, RH_data, mask_nan=mask_nan)
    return prop_asymmetry
