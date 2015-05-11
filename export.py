"""
Export derived measures spreadsheet
(useful for upload to the data exploration tool)
"""
import sys

import numpy as np

from ping.data import PINGData
from research.grouping import group_and_execute

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


def do_usage(args):
    print("\nUsage for %s:" % args[0])
    print("\t%s [prefixes] [groupings]" % args[0])
    print("\n\tExports data based on filtering and grouping;")
    print("\tsetting filters and groups via command-line NYI.")
    print("\n\tprefixes: (optional) comma-delimited list of prefixes")
    print("\t\tto filter computations/results to.")
    print("\tgropuings: (optional) comma-delimited list of ways to")
    print("\t\tsplit the data into groups. A CSV file will be output")
    print("\t\tfor every combination of group unique values.")


if __name__ == '__main__':

    # Filtering / grouping defaults
    filter_args = {
        'prefixes': PINGData.IMAGING_PREFIX,
        'groupings': ['FDH_23_Handedness_Prtcpnt'],
        'limits': {}}
    #    'MRI_cort_area_ctx_frontalpole_AI':
    #        lambda vals: select_lowest_values(-vals)}

    # Parse args
    n_args = len(sys.argv)
    if n_args >= 2:
        filter_args['prefixes'] = sys.argv[1].split(',')
    if n_args >= 3:
        filter_args['groupings'] = sys.argv[2].split(',')
    if n_args >= 4:
        do_usage(sys.argv)

    group_and_execute(fn=export_data, **filter_args)

