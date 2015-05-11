"""
"""
import copy
from functools import partial

import numpy as np

from .data import get_derived_data, get_all_data


def select_lowest_values(vals, pct=0.25):
    """One possible 'limit' filter."""
    selected_idx = np.argsort(vals)[:int(np.floor(pct * len(vals)))]
    out_idx = np.zeros((len(vals),), dtype=bool)
    out_idx[selected_idx] = True
    return out_idx


def group_and_execute(fn, all_data=None, prefixes=None, groupings=None, limits=None):
    """Filters and/or groups data, then runs the given function."""

    # Massage inputs because Python sucks at default args
    #   for lists and dicts.
    if prefixes is None:
        prefixes = []
    if groupings is None:
        groupings = []
    if limits is None:
        limits = dict()

    # Convert prefixes to filter functions
    filter_fns = dict([(p, partial(lambda k, v, p: k.startswith(p), p=p))
                       for p in prefixes])


    if not groupings and not limits:
        print("Case 1: Apply %s filters, return the result." % len(filter_fns))
        data = get_all_data(all_data=all_data, filter_fns=filter_fns)
        return fn(data)
        # data.export(out_file=EXPORTED_PING_SPREADSHEET)

    elif not groupings:
        filter_fns.update(limits)

        print("Case 2: Apply %s filters, return the result." % len(filter_fns))
        data = get_all_data(all_data=all_data, filter_fns=filter_fns)
        return fn(data, limits=limits)
        # filter_and_export(data, limits=limits)

    else:
        # Case 3: do the filtering, apply the grouping, and sub-filter
        #   the group data.
        print("Case 3: Apply %s filters, %s groupings, %s filters, return the result." % (
            len(filter_fns), len(groupings), len(limits)))

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
                fn(data, limits=cur_limits, group={group_key: group_val})
