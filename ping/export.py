"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import os

import csv
import matplotlib.pyplot as plt
import numpy as np

from .access import load_PING_data, get_twohemi_keys
from .utils import do_and_plot_regression

EXPORTED_PING_SPREADSHEET = 'csv/PING.csv'


def compute_all_totals(prefix):
    """ Loop over all properties to show asymmetry."""
    data = load_PING_data(scrub_fields=False)

    export_data = dict((('SubjID', data['SubjID'],),))

    # Process & plot the data.
    for prop_name in get_twohemi_keys(prefix, data.keys()):
        lh_prop_name = prop_name.replace('_rh_', '_lh_').replace('_Right_', '_Left_').replace('_R_', '_L_')
        dest_prop_name = prop_name.replace('_rh_', '_').replace('_Right_', '_').replace('_R_', '_')
        dest_prop_name += '_LH_PLUS_RH'
        export_data[dest_prop_name] = data[prop_name] + data[lh_prop_name]

    return export_data


def combine_genetic_data(export_data, gene_file):
    # Merge the two together.
    with open(gene_file, 'rb') as fp:
        gene_data = np.recfromcsv(fp)

    all_idx = []
    for subj_id in export_data['SubjID']:
        x_subj_id = '"%s"' % subj_id
        if x_subj_id in  gene_data['subjid']:
            all_idx.append(gene_data['subjid'].tolist().index(x_subj_id))
        else:
            all_idx.append(np.nan)

    for key in gene_data.dtype.names:
        if key == 'subjid':
            continue
        vals = [gene_data[idx] if np.logical_not(np.isnan(idx)) else np.nan
                for idx in all_idx]
        export_data[key] = vals

    return export_data


def get_all_derived_data(prefix=None, force=True):
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


def export_all_data(export_data, out_file):

    keys = export_data.keys()

    # Now export as csv
    with open(out_file, 'wb') as fp:
        w = csv.writer(fp)
        w.writerow(keys)
        for row_idx in range(len(export_data['SubjID'])):
            row = []
            for key in keys:
                row.append(export_data[key][row_idx])
            w.writerow(row)
