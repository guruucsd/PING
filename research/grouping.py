"""
"""
import copy
import sys
from functools import partial

import numpy as np

from .data import get_derived_data, get_all_data
from ping.data import PINGData


def select_lowest_values(vals, pct=0.25):
    """One possible 'limit' filter."""
    selected_idx = np.argsort(vals)[:int(np.floor(pct * len(vals)))]
    out_idx = np.zeros((len(vals),), dtype=bool)
    out_idx[selected_idx] = True
    return out_idx


def get_groupings(data, grouping_keys):
    """ Group data by the unique values of `grouping_key`"""

    if not isinstance(grouping_keys, list):
        grouping_keys = [grouping_keys]

    # Group
    grouping_index = group_names = []
    for gpn in grouping_keys:
        grouping_data = np.asarray(data[gpn].tolist())
        prop_groups = (set(np.unique(grouping_data)) -
                       set(['Not yet established', 'nan']))
        prop_groups = np.asarray(list(prop_groups))
        n_prop_group = len(prop_groups)

        if len(group_names) == 0:
            group_names = prop_groups.tolist()
            grouping_index = [grouping_data == pg for pg in prop_groups]
        else:
            group_names = ['%s_%s' % (g1, pg)
                           for g1 in group_names
                           for pg in prop_groups]
            grouping_index = [np.logical_and(idx1, grouping_data == pg)
                              for idx1 in grouping_index
                              for pg in prop_groups]

    return group_names, grouping_index


def select_group_data(data, idx, keys=None, remove_nan=False):
    if keys is None:
        keys = data.keys()

    out_data = dict()
    for key in keys:
        out_data[key] = data[key][idx]

    if remove_nan and len(out_data) > 0:
        good_idx = np.ones(next(iter(out_data.values())).shape, dtype=bool)
        for key, val in out_data.items():
            try:
                good_idx = np.logical_and(good_idx, ~np.isnan(val))
            except TypeError:
                print("Eliminating field %s" % key)
        for key in out_data:
            out_data[key] = out_data[key][good_idx]

    return out_data


def group_and_compare(fn,
                      grouping_keys,
                      data=None,
                      prefixes=None,
                      limits=None,
                      symmetric=True,
                      **kwargs):
    """ Loop over all properties to show asymmetry."""
    # Massage inputs because Python sucks at default args
    #   for lists and dicts.
    if data is None:

        # Convert prefixes to filter functions
        filter_fns = dict([(p, partial(lambda k, v, p: k.startswith(p), p=p))
                           for p in prefixes])
        data = get_all_data(filter_fns=filter_fns).data_dict

    if prefixes is None:
        prefixes = []
    if limits is None:
        limits = dict()

    # Process & plot the data.
    argout = []

    group_names, grouping_index = get_groupings(data, grouping_keys)

    # Collect the data
    group_data = dict()
    for gi, group_name in enumerate(group_names):

        # Index the current group
        if group_name == 'all':
            idx = np.ones((len(data.values()[0]),), dtype=bool)
        else:
            idx = grouping_index[gi]

        # Select data within the group
        group_data[group_name] = select_group_data(data, idx, 
                                                   keys=kwargs.get('keys'),
                                                   remove_nan=kwargs.get('remove_nan', False))

    def group_info(gi1, gi2, curi, group_names):
        n_groups = len(group_names)
        n_total = n_groups * (n_groups - 1) / 2 if symmetric else n_groups**2
        group_info = {
            'gi1': gi1,
            'gi2': gi2,
            'total_compares': n_total,
            'curi': curi,
            'n_groups': n_groups,
            'group_names': group_names
        }
        return group_info

    results = []
    if symmetric:
        for gi1, gn1 in enumerate(group_names):
            for gi2, gn2 in enumerate(group_names[(gi1 + 1):]):
                print(gi1, gi2)
                results.append(fn(group_data[gn1], group_data[gn2], 
                                  group_info(gi1=gi1, gi2=gi1+1+gi2, curi=len(results), 
                                             group_names=group_names),
                                  **kwargs))

    else:
        for gi1, gn1 in enumerate(group_names):
            for gi2, gn2 in enumerate(group_names):
                results.append(fn(group_data[gn1], group_data[gn2], 
                                  group_info(gi1=gi1, gi2=gi2, curi=len(results), 
                                             group_names=group_names),
                                  **kwargs))

    return group_names, results


def group_and_execute(fn, all_data=None, prefixes=None, groupings=None, limits=None, verbose=0, remove_nan=False, **kwargs):
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
    if remove_nan:
        def no_nan_combo(k, v, fn):
            """Return that it's not nan, and some extra function"""
            if k == 'SubjID':
                return True
            try:
                is_nan = np.isnan(v)
                fn_result = fn(k, v)
            except TypeError:
                return False  # mark non-numeric columns to remove.

            if np.all(is_nan):
                return False # remove columns full of NaN to remove.
            elif isinstance(fn_result, bool):
                if not fn_result:
                    return fn_result
                else:
                    return ~is_nan
            else:
                return np.logical_and(fn_result, ~is_nan)  # remove nan elements

        for key, filter_fn in filter_fns.items():
            filter_fns[key] = partial(no_nan_combo, fn=filter_fn)

    if not groupings and not limits:
        print("Case 1: Apply %s filters, return the result." % len(filter_fns))
        data = get_all_data(all_data=all_data, filter_fns=filter_fns, verbose=verbose)
        return fn(data, **kwargs)

    elif not groupings:
        filter_fns.update(limits)

        print("Case 2: Apply %s filters, return the result." % len(filter_fns))
        data = get_all_data(all_data=all_data, filter_fns=filter_fns, verbose=verbose)
        return fn(data, limits=limits, **kwargs)
        # filter_and_export(data, limits=limits)

    else:
        # Case 3: do the filtering, apply the grouping, and sub-filter
        #   the group data.
        print("Case 3: Apply %s filters, %s groupings, %s filters, return the result." % (
            len(filter_fns), len(groupings), len(limits)))

        # Get the data (first pass), for filtering.
        print("Computing data for unfiltered groups...")
        data = get_all_data(filter_fns=filter_fns, verbose=verbose)

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
                group_data = get_all_data(all_data=group_data, filter_fns=filter_fns, verbose=verbose)

                # Now export
                fn(data, limits=cur_limits, group={group_key: group_val}, **kwargs)


def parse_filter_args(args, filter_defaults=None):

    # Filtering / grouping defaults
    filter_args = {
        'prefixes': PINGData.IMAGING_PREFIX,
        'groupings': [],
        'limits': {}}
    #    'MRI_cort_area_ctx_frontalpole_AI':
    #        lambda vals: select_lowest_values(-vals)}

    # Update defaults by those passed in
    if filter_defaults is not None:
        filter_args.update(filter_defaults)

    # Parse args
    n_args = len(args)
    if n_args >= 1:
        filter_args['prefixes'] = args[0].split(',')
    if n_args >= 2:
        filter_args['groupings'] = args[1].split(',')
    return filter_args


def chunk_string(str, max_len=80):
    """breaks a string into substrings below 
    """
    chunks = ['']
    for w in str.split(' '):
        extend_chunk = ' '.join([chunks[-1], w])
        if len(extend_chunk) <= max_len:
            chunks[-1] = extend_chunk
        else:
            chunks.append(w.strip())

    return chunks


def do_usage_grouping(exec_name, description="", args=None, optargs=None, error_msg=None):
    """Print out usage."""

    if error_msg:
        print("*** ERROR: %s" % error_msg)
    print("\n%s %s[prefixes] [groupings]" % (exec_name, " ".join(args or [])))

    if description:
        chunks = chunk_string(description, max_len=60)
        print("\t%s" % chunks[0])
        for chunk in chunks[1:]:
            print("\t\t%s" % chunk)

    if args is not None:
        print("")
        for arg_name, arg_desc in args.items():
            chunks = chunk_string(arg_desc, max_len=60)
            print("%s: %s" % (arg_name, chunks[0]))
            for chunk in chunks[1:]:
                print("\t%s" % chunk)

    print("\nprefixes: (optional) comma-delimited list of prefixes")
    print("\tto filter computations/results to.")
    print("groupings: (optional) comma-delimited list of ways to")
    print("\tsplit the data into groups. A CSV file will be output")
    print("\tfor every combination of group unique values.")

    if optargs is not None:
        for arg_name, arg_desc in optargs.items():
            chunks = chunk_string(arg_desc, max_len=55)
            print("%s: (optional) %s" % (arg_name, chunks[0]))
            for chunk in chunks[1:]:
                print("\t%s" % chunk)

