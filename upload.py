"""
Export all measures to a spreadsheet, then
upload to the PING server.
"""
from functools import partial

from export import EXPORTED_PING_SPREADSHEET
from ping.apps import PINGSession
from research.data import get_all_data


def do_upload(*args):
    if len(args) > 0:
        # Someone told us about an existing csv
        csv_file = args[0]

    else:
        # We will build a csv on the fly.
        csv_file = EXPORTED_PING_SPREADSHEET

        # Limit to imaging data only
        data = get_all_data()
        data.filter(filter_fns=[partial(lambda k, v, p: k.startswith(p), p=p)
                                for p in data.IMAGING_PREFIX])

        data.export(out_file=csv_file)

    # Upload the new spreadsheet.
    sess = PINGSession()
    sess.login()
    sess.upload_user_spreadsheet(csv_file)


if __name__ == '__main__':
    import sys
    do_upload(*sys.argv[1:])
