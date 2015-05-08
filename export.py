"""
Export
"""
from data import get_derived_data


if __name__ == '__main__':
    prefixes = ['MRI_cort_area', 'MRI_cort_thick', 'MRI_subcort_vol', 'DTI_fiber_vol']
    export_data = get_derived_data(prefix=prefixes, force=True)
    export_all_data(export_data, out_file=EXPORTED_PING_SPREADSHEET)
