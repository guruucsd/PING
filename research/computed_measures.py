"""
"""
import collections
import copy
import os

import numpy as np
import pandas

from .asymmetry import (get_asymmetry_index, get_ai_key,
                        is_ai_key)
from .multivariate import AsymmetryPCA
from ping.data import PINGData, get_lh_key, get_nonhemi_key


def compute_all_totals(filtered_data):
    """ Loop over all properties to show asymmetry."""

    # Process & plot the data.
    out_data = PINGData(dict(SubjID=filtered_data.data_dict['SubjID']))
    for key in filtered_data.get_twohemi_keys():
        lh_key = get_lh_key(key)
        dest_key = get_nonhemi_key(key)
        dest_key += '_LH_PLUS_RH'
        out_data.data_dict[dest_key] = filtered_data.data_dict[key] + filtered_data.data_dict[lh_key]

    return out_data


def compute_all_asymmetries(filtered_data):
    """ Loop over all properties to show asymmetry."""
    # Process & plot the data.
    out_data = PINGData(dict(SubjID=filtered_data.data_dict['SubjID']))
    for key in filtered_data.get_twohemi_keys():
        dest_key = get_ai_key(key)
        out_data.data_dict[dest_key] = get_asymmetry_index(filtered_data.data_dict,
                                                           key, mask_nan=False)

    # Add one property for total asymmetry
    n_subj = out_data.get_num_subjects()
    for p in filtered_data.IMAGING_PREFIX:
        # Compute total asymmetry per subject
        total_asymmetry = np.zeros((n_subj,))
        good_keys = filter(lambda k: (k.startswith(p) and 
                                      is_ai_key(k)),
                           out_data.data_dict.keys())
        for key in good_keys:
            values = out_data.data_dict[key].copy()
            values[np.isnan(values)] = 0. # summing, so is ok.
            total_asymmetry += out_data.data_dict[key]**2

        total_ai_key = get_ai_key(p + '_TOTAL')
        out_data.data_dict[total_ai_key] = np.sqrt(total_asymmetry)

    return out_data


def compute_component_loadings(filtered_data, verbose=0):
    asymmetry_data = compute_all_asymmetries(filtered_data)

    pca = AsymmetryPCA(whiten=True)
    pca.fit(asymmetry_data, verbose=verbose)

    if verbose >= 1:
        pca.report_asymmetry_loadings()
        pca.report_behavior_correlations()
        pca.report_background_correlations()

    pc_projections = pca.get_projections()
    keys = ['PC%d' % pi for pi in range(pc_projections.shape[0])]

    pc_dict = dict(zip(keys, pc_projections))
    pc_dict['SubjID'] = pca.subj_ids

    return PINGData(pc_dict)


def get_derived_data(filtered_data, tag=None, verbose=0):
    print("Computing derived data...")

    data = compute_all_asymmetries(filtered_data)
    data.merge(compute_all_totals(filtered_data))

    # Add principle components as a whole,
    #   then for each prefix seperately.
    cl = compute_component_loadings(filtered_data, verbose=verbose)
    data.merge(cl, tag=tag)
    return data


def get_all_data(all_data=None, filter_fns=None, verbose=0):
    if all_data is None:
        all_data = PINGData()

    if filter_fns is None:
        filter_fns = dict()

    # Now get derived data for each filter separately,
    # and all filters together.
    for key in ['ALL'] + list(filter_fns.keys()):
        filtered_data = copy.deepcopy(all_data)
        if key == 'ALL':
            filtered_data.filter(list(filter_fns.values()))
        else:
            filtered_data.filter(filter_fns[key])
        computed_data = get_derived_data(filtered_data, tag='%s_%%s' % key,
                                         verbose=verbose)
        all_data.merge(computed_data)

    return all_data
