"""
"""
import collections
import copy
import inspect
import os

import numpy as np
import pandas
from six import string_types

from ..asymmetry import (get_asymmetry_index, get_ai_key,
                         is_ai_key)
from ..multivariate import AsymmetryPCA
from ...ping.data import PINGData, DestrieuxData


def keytype2label(key):
    if not isinstance(key, string_types):
        return [keytype2label(k) for k in key]

    if key.endswith('TOTAL_LH_PLUS_RH'):
        return 'Total'
    elif key.endswith('LH_PLUS_RH'):
        return 'Mean'
    elif key.endswith('TOTAL_AI'):
        return 'Overall asymmetry index'
    elif key.endswith('AI'):
        return 'Asymmetry Index'
    else:
        return key

def compute_all_totals(filtered_data):
    """Computes total surface area / volume."""

    # Process & plot the data.
    out_data = filtered_data.__class__(dict(SubjID=filtered_data.data_dict['SubjID']))
    for rh_key in filtered_data.get_twohemi_keys():
        assert filtered_data.which_hemi(rh_key) == 'rh'
        lh_key = filtered_data.get_lh_key_from_rh_key(rh_key)
        dest_key = filtered_data.get_nonhemi_key(rh_key)
        dest_key += '_LH_PLUS_RH'
        out_data.data_dict[dest_key] = filtered_data.data_dict[rh_key] + filtered_data.data_dict[lh_key]

    def is_sum_key(key):
        return key.endswith('_LH_PLUS_RH')
    def get_sum_key(key):
        return '%s_LH_PLUS_RH' % key

    # Add one property for total asymmetry
    n_subj = out_data.get_num_subjects()
    for p in filtered_data.IMAGING_PREFIX:
        good_keys = filter(lambda k: (k.startswith(p) and
                                      is_sum_key(k)),
                           out_data.data_dict.keys())
        good_keys = list(good_keys)
        if len(good_keys) == 0:
            continue

        # Compute total LH+RH, per subject
        total_area = np.zeros((n_subj,))
        total_measures = np.zeros((n_subj,))
        for key in good_keys:
            values = out_data.data_dict[key]
            good_idx = np.logical_not(np.isnan(values))
            total_area[good_idx] += values[good_idx]**2
            total_measures[good_idx] += 1
        total_sum_key = get_sum_key(p + '_TOTAL')
        out_data.data_dict[total_sum_key] = np.sqrt(total_area / total_measures)
    return out_data


def compute_all_volumes(filtered_data):
    """ Loop over all properties to show asymmetry."""

    # Process & plot the data.
    out_data = filtered_data.__class__(dict(SubjID=filtered_data.data_dict['SubjID']))
    for rh_key_area in [k for k in filtered_data.get_twohemi_keys()
                          if k.startswith('MRI_cort_area.ctx')]:
        assert filtered_data.which_hemi(rh_key_area) == 'rh'
        try:
            rh_key_thick = rh_key_area.replace('MRI_cort_area.ctx', 'MRI_cort_thick.ctx')
            rh_key_vol = rh_key_area.replace('MRI_cort_area.ctx', 'MRI_cort_vol.ctx')
            lh_key_area = filtered_data.get_lh_key_from_rh_key(rh_key_area)
            lh_key_thick = filtered_data.get_lh_key_from_rh_key(rh_key_thick)
            lh_key_vol = filtered_data.get_lh_key_from_rh_key(rh_key_vol)
            dest_key = filtered_data.get_nonhemi_key(rh_key_vol)
            out_data.data_dict[dest_key] = filtered_data.data_dict[rh_key_vol] + filtered_data.data_dict[lh_key_vol]
        except Exception as e:
            print(e)
            continue

    return out_data


def compute_all_asymmetries(filtered_data):
    """ Loop over all properties to show asymmetry."""

    # Asymmetry index
    out_data = filtered_data.__class__(dict(SubjID=filtered_data.data_dict['SubjID']))

    for key in filtered_data.get_twohemi_keys():
        dest_key = get_ai_key(key, filtered_data)
        out_data.data_dict[dest_key] = get_asymmetry_index(filtered_data,
                                                           key, mask_nan=False)

    # Add one property for total asymmetry
    n_subj = out_data.get_num_subjects()
    for p in filtered_data.IMAGING_PREFIX:
        good_keys = filter(lambda k: (k.startswith(p) and
                                      is_ai_key(k)),
                           out_data.data_dict.keys())
        good_keys = list(good_keys)
        if len(good_keys) == 0:
            continue

        # Compute total asymmetry per subject
        total_asymmetry = np.zeros((n_subj,))
        total_measures = np.zeros((n_subj,))
        for key in good_keys:
            values = out_data.data_dict[key]
            good_idx = np.logical_not(np.isnan(values))
            total_asymmetry[good_idx] += values[good_idx]**2
            total_measures[good_idx] += 1
        total_ai_key = get_ai_key(p + '%s', filtered_data) % '_TOTAL'
        out_data.data_dict[total_ai_key] = np.sqrt(total_asymmetry / total_measures)

    return out_data


def compute_component_loadings(filtered_data, verbose=0):
    # Volume before asymmetry, so asymmetry can be computed.
    volume_data = compute_all_volumes(filtered_data)
    asymmetry_data = compute_all_asymmetries(filtered_data)

    pca = AsymmetryPCA(whiten=False)
    pca.fit(asymmetry_data, verbose=verbose)

    if verbose >= 1:
        pca.report_asymmetry_loadings()
        pca.report_behavior_correlations()
        pca.report_background_correlations()

    pc_projections = pca.get_projections()
    keys = ['PC%d' % pi for pi in range(pc_projections.shape[0])]

    pc_dict = dict(zip(keys, pc_projections))
    pc_dict['SubjID'] = pca.subj_ids

    return filtered_data.__class__(pc_dict)


def get_derived_data(filtered_data, tag=None, verbose=0):
    print("Computing derived data...")

    data = compute_all_asymmetries(filtered_data)
    data.merge(compute_all_totals(filtered_data))

    # Add principle components as a whole,
    #   then for each prefix seperately.
    cl = compute_component_loadings(filtered_data, verbose=verbose)
    data.merge(cl, tag=tag)
    return data


known_data = dict(desikan=dict(klass=PINGData),
                  destrieux=dict(klass=DestrieuxData))

def get_all_data(all_data='desikan', filter_fns=None, verbose=0,
                 username=None, passwd=None, data_dir=None):
    kwargs = dict(username=username, passwd=passwd)
    if data_dir is not None:
        kwargs['data_dir'] = data_dir

    if inspect.isclass(all_data):
        data_klass = all_data
        all_data = data_klass(**kwargs)
    elif isinstance(all_data, string_types):
        dataset_name = all_data
        global known_data
        if dataset_name not in known_data:
            raise ValueError('Unknown dataset: %s' % dataset_name)
        else:
            if not known_data[dataset_name].get('data'):
                data_klass = known_data[dataset_name]['klass']
                known_data[dataset_name]['data'] = data_klass(**kwargs)
            all_data = known_data[dataset_name]['data']
    else:
        raise NotImplementedError()  # assume it's data.

    if filter_fns is None:
        filter_fns = dict()

    # Now get derived data for each filter separately,
    # and all filters together.
    out_data = copy.deepcopy(all_data)  # PINGData(dict(SubjID=all_data.data_dict['SubjID']))
    for key in ['ALL'] + list(filter_fns.keys()):
        filtered_data = copy.deepcopy(all_data)
        if key == 'ALL':
            filters = list(filter_fns.values())
        else:
            filters = filter_fns[key]

        filtered_data.filter(filters)
        computed_data = get_derived_data(filtered_data, tag='%s_%%s' % key,
                                         verbose=verbose)
        computed_data.filter(filters)
        out_data.merge(computed_data)

    # Now filter out any subjects with all nan measures.
    out_data.purge_empty_subjects()

    return out_data
