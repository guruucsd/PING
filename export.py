"""
Export derived measures spreadsheet
(useful for upload to the data exploration tool)
"""
import copy
from functools import partial

import numpy as np

from ping.data import PINGData
from research.computed_measures import get_derived_data, get_all_data

EXPORTED_PING_SPREADSHEET = 'csv/PING_userdefined.csv'


def load_user_spreadsheet():
    print("Loading derived data...")
    data = pandas.read_csv(EXPORTED_PING_SPREADSHEET, low_memory=False)
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


def filter_and_export(data, limits, group=None, group_val=None):
    """"""
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

    # Filter the data
    group_data = copy.deepcopy(data)

    group_data.filter(list(limits.values()))
    
    # Export
    group_data.export(out_file=cur_csv)

    return cur_csv


if __name__ == '__main__':
    prefixes = PINGData.IMAGING_PREFIX
    groupings = ['FDH_23_Handedness_Prtcpnt']
    limits = {}
    #    'MRI_cort_area_ctx_frontalpole_AI':
    #        lambda vals: select_lowest_values(-vals)}

    filter_fns = dict([(p, partial(lambda k, v, p: k.startswith(p), p=p))
                       for p in prefixes])

    if not groupings and not limits:
        data = get_all_data(filter_fns)
        data.export(out_file=EXPORTED_PING_SPREADSHEET)

    elif not groupings:
        data = get_all_data(filter_fns=filter_fns)
        filter_and_export(data, limits=limits)

    else:
        # Get the data (first pass), for filtering.
        print("Computing data for unfiltered groups...")
        data = get_all_data(filter_fns=filter_fns)

        for group_key in groupings:
            group_vals = np.asarray(data.data_dict[group_key])

            for group_val in np.unique([str(v) for v in data.data_dict[group_key]]):
                # Set filters
                cur_limits = copy.deepcopy(limits)
                cur_limits[group_key] = partial(lambda k, v, gk, gv: k != gk or v == gv,
                                                gk=group_key, gv=group_val)

                # Filter the data
                group_data = copy.deepcopy(data)
                group_data.filter(list(cur_limits.values()))
                if group_data.get_num_subjects() == 0:
                    print("Skipping empty group %s..." % group_val)
                    continue

                # Recompute derived data based on group.
                print("Recomputing data for group %s..." % group_val)
                group_data = get_all_data(all_data=group_data, filter_fns=filter_fns)

                # Now export
                filter_and_export(data, limits=cur_limits,
                                  group=group_key, group_val=group_val)
