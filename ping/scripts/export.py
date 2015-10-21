"""
Export derived measures spreadsheet
(useful for upload to the data exploration tool)
"""
import sys

from ..ping.apps import PINGSession
from ..research.grouping import (do_usage_grouping)

EXPORTED_PING_SPREADSHEET = 'data/PING_userdefined.csv'


def parse_filter_args(args, filter_defaults=None):
    warnings.warn('Still using PINGData; must migrate!', DeprecationWarning)

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


def group_and_execute(fn, all_data='desikan', prefixes=None, groupings=None,
                      limits=None, verbose=0, remove_nan=False, **kwargs):
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


def select_lowest_values(vals, pct=0.25):
    """One possible 'limit' filter."""
    selected_idx = np.argsort(vals)[:int(np.floor(pct * len(vals)))]
    out_idx = np.zeros((len(vals),), dtype=bool)
    out_idx[selected_idx] = True
    return out_idx


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


def do_upload(csv_file=EXPORTED_PING_SPREADSHEET,
              atlas='desikan', data_dir='data',
              username=None, passwd=None, force=False):

    if force or not os.path.exists(csv_file):
        # Limit to imaging data only
        data = get_all_data(atlas, data_dir=data_dir,
                            username=username, passwd=passwd)
        # data.filter(filter_fns=[partial(lambda k, v, p: k.startswith(p), p=p)
        #                         for p in data.IMAGING_PREFIX])
        data.export(out_file=csv_file)

    # Upload the new spreadsheet.
    sess = PINGSession()
    sess.login()
    sess.upload_user_spreadsheet(csv_file)


def get_csv_filename(base_filename=EXPORTED_PING_SPREADSHEET, limits=None, group=None):
    """"""
    # Compute a text tag
    tags = []
    if group:
        tags.append('__'.join(['%s_eq_%s' % (key, val.replace(' ', ''))
                               for key, val in group.items()]))
    if limits:
        tags.append('__'.join(limits.keys()))
    tag = '_limit_'.join(tags)

    # Compute a CSV
    cur_csv = base_filename
    if tag:
        cur_csv = cur_csv.replace('.csv', '_%s.csv' % tag)
    print("Dumping data to %s" % (cur_csv))

    return cur_csv


def export_data(data, **kwargs):
    cur_csv = get_csv_filename(**kwargs)
    data.export(out_file=cur_csv)
    return cur_csv


def do_export(prefixes, groupings=None,
              atlas='desikan', data_dir='data',
              username=None, passwd=None,
              force=False):
    filter_args = parse_filter_args(args)
    group_and_execute(fn=export_data, verbose=1, **filter_args)


if __name__ == '__main__':
    parser = ResearchArgParser(description="Upload data spreadsheet to the PING portal.",
                               common_args=['prefixes', 'atlas', 'data-dir',
                                            'username', 'passwd', 'force'])
    parser.add_argument('--groupings', nargs='?' default=None,
                        choices=['Gender', 'FDH_23_Handedness_Prtcpnt'])
    parser.add_argument('--csv_file', nargs='?', default=EXPORTED_PING_SPREADSHEET)

    args = parser.parse_args()
    do_export(**vars(args))

