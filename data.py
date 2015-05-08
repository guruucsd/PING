from ping.access import get_lh_prop_name, get_nonhemi_prop_name
from ping.asymmetry import compute_all_asymmetries
from ping.export import export_all_data

EXPORTED_PING_SPREADSHEET = 'csv/PING.csv'


def compute_all_totals(prefix):
    """ Loop over all properties to show asymmetry."""
    data = load_PING_data(scrub_fields=False)

    export_data = dict((('SubjID', data['SubjID'],),))

    # Process & plot the data.
    for prop_name in get_twohemi_keys(data.keys(), prefix=prefix):
        lh_prop_name = get_lh_prop_name(prop_name)
        dest_prop_name = get_nonhemi_prop_name(prop_name)
        dest_prop_name += '_LH_PLUS_RH'
        export_data[dest_prop_name] = data[prop_name] + data[lh_prop_name]

    return export_data


def get_derived_data(prefix=None, csv_files=[], force=True):
    prefix = prefix or []
    data = None
    if os.path.exists(EXPORTED_PING_SPREADSHEET) and not force:
        print("Loading derived data...")
        data = pandas.read_csv(EXPORTED_PING_SPREADSHEET)
        for p in prefix:
            if not np.any([key.startswith(p) for key in data.keys()]):
                data = None
                break
        if data is not None:
            new_data = dict()
            for key in data.keys():
                new_data[key.replace('.', '_')] = data[key]
            data = new_data

    if data is None:
        print "Computing derived data..."
        data = compute_all_asymmetries(prefix=prefix)
        data.update(compute_all_totals(prefix=prefix))
        data = combine_genetic_data(data, 'csv/frontalpole_genes.csv')

    return data


def get_all_data(prefix):
    all_data = get_derived_data(prefix=prefix)
    ping_data = copy.deepcopy(load_PING_data())

    for key in copy.deepcopy(ping_data.keys()):
        if not np.any([key.startswith(p) for p in prefix]):
            del ping_data[key]
    all_data.update(ping_data)
    return all_data
