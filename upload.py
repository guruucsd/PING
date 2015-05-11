"""
Export all measures to a spreadsheet, then
upload to the PING server.
"""
import sys
from functools import partial

from export import EXPORTED_PING_SPREADSHEET
from ping.apps import PINGSession
from ping.data import PINGData
from research.computed_measures import get_all_data


if len(sys.argv) > 1:
    # Someone told us about an existing csv
    csv_file = sys.argv[1]

else:
    # We will build a csv on the fly.
    csv_file = EXPORTED_PING_SPREADSHEET
    filter_fns = [partial(lambda k, v: k.startswith(p), p=p)
                  for p in PINGData.IMAGING_PREFIX]

    data = get_all_data(filter_fns=filter_fns)
    data.export(out_file=csv_file)

    for key in sorted(data.data_dict.keys()):
        print(key)

# Upload the new spreadsheet.
sess = PINGSession()
sess.login()
sess.upload_user_spreadsheet(csv_file)
