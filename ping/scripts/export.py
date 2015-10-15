"""
Export derived measures spreadsheet
(useful for upload to the data exploration tool)
"""
import sys

from ..ping.apps import PINGSession
from ..research.grouping import (do_usage_grouping, group_and_execute,
                                 parse_filter_args)

EXPORTED_PING_SPREADSHEET = 'data/PING_userdefined.csv'


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


def do_usage(args, error_msg=None):
    do_usage_grouping(__file__, error_msg=error_msg,
                      description=("Exports data based on filtering and grouping; "
                                   "setting filters and groups via command-line NYI."))


def do_export(*args):
    if len(args) > 3:
        do_usage(args, error_msg="Too many arguments.")
        return

    filter_args = parse_filter_args(args)
    group_and_execute(fn=export_data, verbose=1, **filter_args)


if __name__ == '__main__':
    do_export(*sys.argv[1:])
