"""
Export derived measures spreadsheet
(useful for upload to the data exploration tool)
"""
import sys

from research.grouping import (do_usage, group_and_execute,
                               parse_filter_args, UsageError)

EXPORTED_PING_SPREADSHEET = 'csv/PING_userdefined.csv'


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


if __name__ == '__main__':
    try:
        filter_args = parse_filter_args(sys.argv[1:])
    except UsageError:
        do_usage(sys.argv[0], ("Exports data based on filtering and grouping; "
                               "setting filters and groups via command-line NYI."))
    else:
        group_and_execute(fn=export_data, verbose=1, **filter_args)

