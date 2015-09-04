import collections
import copy
import itertools
import warnings

import numpy as np


def filter_data(data_dict, filter_fns, tag=None):
    """Filter functions must take a key and ndarray of values.
    They can return a single boolean (False) to accept/reject the values,
    or an ndarray to select the indices.

    Filters using 'and' logic.
    """
    if filter_fns is None or len(data_dict) == 0:
        return data_dict
    elif not isinstance(filter_fns, dict):
        filter_fns = dict(zip(data_dict.keys(), itertools.cycle([filter_fns])))

    # Limit based on group value and filtering function
    col_len = len(next(iter(data_dict.values())))
    selected_idx = np.ones((col_len,), dtype=bool)
    selected_keys = list(data_dict.keys())
    for filter_key, filter_fn in filter_fns.items():


        if filter_key not in data_dict:
            warnings.warn('%s has a filter, but not found in data; skipping.' % filter_key)
            continue

        # For each key, can have multiple filters
        if not isinstance(filter_fn, collections.Iterable):
            filter_fn = [filter_fn]

        for fn in filter_fn:
            # If we got here, we have a valid key and a valid fn
            filter_vals = data_dict[filter_key]
            filter_output = fn(filter_key, filter_vals[selected_idx])
            if isinstance(filter_output, bool) or isinstance(filter_output, np.bool_):
                if not filter_output:
                    selected_keys.remove(filter_key)
                break  # accepted as whole, skip below logic for efficiency.
            new_idx = np.zeros(selected_idx.shape)
            new_idx[selected_idx] = filter_output
            selected_idx = np.logical_and(selected_idx, new_idx)
            assert np.any(selected_idx), 'All items were filtered out.'

    # Pull out filtered values and 'tag' the key
    # (this makes documenting results much easier)
    filtered_data = dict()
    for key in selected_keys:
        if key in ['SubjID'] or tag is None:
            tagged_key = key
        else:
            tagged_key = '%s_%s' % (key, tag)

        filtered_data[tagged_key] = np.asarray(data_dict[key])[selected_idx]

    return filtered_data


def merge_by_key(dict1, dict2, merge_key='SubjID', tag=None):

    # Make index mapping from dict1 into dict2
    all_idx = []
    for val in dict1[merge_key]:
        if val in dict2[merge_key]:
            if not isinstance(dict2[merge_key], list):
                dict2[merge_key] = dict2[merge_key].tolist()
            all_idx.append(dict2[merge_key].index(val))
        else:
            all_idx.append(np.nan)
    all_idx = np.asarray(all_idx)

    # Now reorder dict2 values to dict1
    out_dict = copy.deepcopy(dict1)
    for key, vals in dict2.items():
        # Don't reassign primary key
        if key == merge_key:
            continue

        # Reorder dict2's values into dict1's order
        if np.any(np.isnan(all_idx)):
            reordered_vals = np.asarray([vals[idx] if not np.isnan(idx) else np.nan
                                         for idx in all_idx])
        else:
            reordered_vals = vals[all_idx]

        # Assign, possibly after tagging keys
        if tag is not None:
            key = tag % key
        if key in out_dict:
            pass  # warnings.warn('Overwriting "%s" values from original dict in merge.' % key)
        out_dict[key] = reordered_vals

    return out_dict
