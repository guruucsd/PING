"""
Export derived measures spreadsheet
(useful for upload to the data exploration tool)
"""
import copy

import numpy as np

from ping.access import load_PING_data
from ping.computed_measures import get_derived_data
from ping.export import (export_all_data, merge_by_key)

EXPORTED_PING_SPREADSHEET = 'csv/PING_userdefined.csv'


def load_user_spreadsheet():
    print("Loading derived data...")
    data = pandas.read_csv(EXPORTED_PING_SPREADSHEET)
    new_data = dict()
    for key in data.keys():
        new_data[key.replace('.', '_')] = data[key].as_matrix()
    data = new_data
    return data



def select_lowest_values(vals, pct=0.25):
    selected_idx = np.argsort(vals)[:int(np.floor(pct * len(vals)))]
    out_idx = np.zeros((len(vals),), dtype=bool)
    out_idx[selected_idx] = True
    return out_idx


def filter_data(merged_data, limits=None):
    limits = limits or dict()

    # Limit based on group value and filtering function
    for limit_key, limit_fn in limits.items():
        limit_vals = merged_data[limit_key]
        limit_idx = limit_fn(limit_vals[selected_idx])
        new_idx = np.zeros(selected_idx.shape)
        new_idx[selected_idx] = limit_idx
        selected_idx = np.logical_and(selected_idx, new_idx)

    # Pull out filtered values and 'tag' the key
    # (this makes documenting results much easier)
    filtered_data = dict()
    for key, val in export_data.items():
        if key in ['SubjID']:
            tagged_key = key
        else:
            tagged_key = '%s_%s' % (key, tag)
        filtered_data[tagged_key] = np.asarray(val)[selected_idx]

    return filtered_data


def filter_and_export(cur_export_data, limits, group=None, group_val=None):
    """"""
    # cur_export_data = filter_derived_data(prefix=prefixes,
    #                                       limits=limits)

    # Compute a text tag
    tags = []
    if group:
        tags.append('%s_eq_%s' % (group, group_val.replace(' ', '')))
    if limits:
        tags.append('__'.join(limits.keys()))
    tag = '_limit_'.join(tags)

    # Compute a CSV
    cur_csv = EXPORTED_PING_SPREADSHEET.replace('.csv', '_%s.csv' % tag)
    print("Dumping group %s=%s to %s" % (group, group_val, cur_csv))

    export_all_data(cur_export_data, out_file=cur_csv)


if __name__ == '__main__':
    prefixes = ['MRI_cort_area', 'MRI_cort_thick',
                'MRI_subcort_vol', 'DTI_fiber_vol']
    groupings = ['FDH_23_Handedness_Prtcpnt']
    limits = {}
    #    'MRI_cort_area_ctx_frontalpole_AI':
    #        lambda vals: select_lowest_values(-vals)}

    if not groupings and not limits:
        export_data = get_derived_data(prefix=prefixes)
        export_all_data(export_data, out_file=cur_csv)

    elif not groupings:
        merged_data = merge_by_key(get_derived_data(prefix=prefixes),
                                   load_PING_data())
        filter_and_export(merged_data, limits=limits)

    else:
        merged_data = merge_by_key(get_derived_data(prefix=prefixes),
                                   load_PING_data())
        for group in groupings:

            group_vals = np.asarray(merged_data[group])

            for group_val in np.unique(merged_data[group]):
                cur_limit = copy.copy(limits)
                cur_limit[group] = lambda vals: np.asarray(vals) == group_val

                filter_and_export(merged_data, limits=cur_limit,
                                  group=group, group_val=group_val)
