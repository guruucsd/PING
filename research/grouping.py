"""
"""
import copy
import sys
from functools import partial

import numpy as np

from .data import get_derived_data, get_all_data
from ping.data import PINGData


class UsageError(Exception):
    pass


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


def parse_filter_args(args, filter_defaults=None):

    # Filtering / grouping defaults
    filter_args = {
        'prefixes': PINGData.IMAGING_PREFIX,
        'groupings': ['FDH_23_Handedness_Prtcpnt'],
        'limits': {}}
    #    'MRI_cort_area_ctx_frontalpole_AI':
    #        lambda vals: select_lowest_values(-vals)}

    # Update defaults by those passed in
    if filter_defaults is not None:
        filter_args.update(filter_defaults)

    # Parse args
    n_args = len(args)
    if n_args >= 1:
        filter_args['prefixes'] = sys.argv[1].split(',')
    if n_args >= 2:
        filter_args['groupings'] = sys.argv[2].split(',')
    if n_args >= 3:
        raise UsageError(args)


def do_usage(exec_name, description=""):
    print("\n%s [prefixes] [groupings]" % exec_name)

    if description:
        chunks = ['']
        for w in description.split(' '):
            extend_chunk = ' '.join([chunks[-1], w])
            if len(extend_chunk) <= 60:
                chunks[-1] = extend_chunk
            else:
                chunks.append(w)

        print("\t%s" % chunks[0][1:])
        for chunk in chunks[1:]:
            print("\t\t%s" % chunk)

    print("\nprefixes: (optional) comma-delimited list of prefixes")
    print("\tto filter computations/results to.")
    print("groupings: (optional) comma-delimited list of ways to")
    print("\tsplit the data into groups. A CSV file will be output")
    print("\tfor every combination of group unique values.")
    print("\n")