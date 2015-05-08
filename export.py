"""
Export derived measures spreadsheet
(useful for upload to the data exploration tool)
"""
from ping.computed_measures import get_derived_data
from ping.export import export_all_data

EXPORTED_PING_SPREADSHEET = 'csv/PING_userdata.csv'


def load_user_spreadsheet():
    print("Loading derived data...")
    data = pandas.read_csv(EXPORTED_PING_SPREADSHEET)
    new_data = dict()
    for key in data.keys():
        new_data[key.replace('.', '_')] = data[key].as_matrix()
    data = new_data
    return data

prefixes = ['MRI_cort_area', 'MRI_cort_thick', 'MRI_subcort_vol', 'DTI_fiber_vol']
export_data = get_derived_data(prefix=prefixes)
export_all_data(export_data, out_file=EXPORTED_PING_SPREADSHEET)
