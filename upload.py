"""
Export all measures to a spreadsheet, then
upload to the PING server.
"""
import sys

from export import EXPORTED_PING_SPREADSHEET
from ping.computed_measures import get_derived_data
from ping.apps.upload import PINGUploadSession
from ping.export import export_all_data


# Export the relevant data
if len(sys.argv) > 1:
    csv_file = sys.argv[1]
else:
    csv_file = EXPORTED_PING_SPREADSHEET
    prefixes = ['MRI_cort_area', 'MRI_cort_thick', 'MRI_subcort_vol', 'DTI_fiber_vol']
    export_data = get_derived_data(prefix=prefixes)
    export_all_data(export_data, out_file=csv_file)

    for key in sorted(export_data.keys()):
        print key

# Upload the new spreadsheet.
sess = PINGUploadSession()
sess.login()
sess.upload_spreadsheet(csv_file)
