"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import os

import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas

from ping import col2prop, load_PING_data, get_twohemi_keys, get_asymmetry_index
from utils import do_and_plot_regression


def compare_group_asymmetry(data, prop_name, grouping_prop_name):
    """ Groups data according to grouping_prop_name, computes
    asymmetry index for prop_name, and graphs."""

    group_names, grouping_index = make_groups(data, grouping_prop_name)

    age_data = np.asarray(data['Age_At_IMGExam'].tolist())
    fh = plt.figure(figsize=(18, 10))
    n_subplots = 1 + len(group_names)
    n_rows = int(np.round(0.75 * np.sqrt(n_subplots)))
    n_cols = int(np.ceil(n_subplots / float(n_rows)))
    for gi, group_name in enumerate(['all'] + group_names):
        # Index the current group
        if group_name == 'all':
            idx = np.ones((data.shape[0],), dtype=bool)
        else:
            idx = grouping_index[gi - 1]

        # Select data within the group
        cur_ages = age_data[idx]
        prop_asymmetry = get_asymmetry_index(data, prop_name, mask_nan=True)[idx]

        # Remove bad data
        bad_idx = np.logical_or(np.isnan(cur_ages), np.isnan(prop_asymmetry))
        good_idx = np.logical_not(bad_idx)
        cur_ages = cur_ages[good_idx]
        prop_asymmetry = prop_asymmetry[good_idx]

        # Plot the result
        ax = fh.add_subplot(n_rows, n_cols, gi + 1)
        do_and_plot_regression(cur_ages, prop_asymmetry, ax=ax,
                               xlabel='Age', ylabel='Asymmetry Index (LH - RH)',
                               title='Group: %s\n%s' % (group_name, prop_name))


def compute_all_asymmetries(prefix):
    """ Loop over all properties to show asymmetry."""
    data = load_PING_data(scrub_fields=False)

    export_data = dict((('SubjID', data['SubjID'],),))

    # Process & plot the data.
    for prop_name in get_twohemi_keys(prefix, data.keys()):
        dest_prop_name = prop_name.replace('_rh_', '_').replace('_Right_', '_').replace('_R_', '_')
        dest_prop_name += '_AI'
        export_data[dest_prop_name] = get_asymmetry_index(data, prop_name, mask_nan=False)

    # Add total asymmetry
    n_subj = len(data['SubjID'])
    for p in prefix:
        total_asymmetry = np.zeros((n_subj,))
        for key in filter(lambda k: k.startswith(p) and k.endswith('_AI'), export_data.keys()):
            values = export_data[key].copy()
            values[np.isnan(values)] = 0.
            total_asymmetry += export_data[key]**2
        export_data[p + '_TOTAL_AI'] = np.sqrt(total_asymmetry)
        assert len(export_data[p + '_TOTAL_AI']) == n_subj
    return export_data


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



if __name__ == '__main__':
    prefixes = ['MRI_cort_area', 'MRI_cort_thick', 'MRI_subcort_vol', 'DTI_fiber_vol']
    export_data = compute_all_asymmetries(prefix=prefixes)
    export_data.update(compute_all_totals(prefix=prefixes))
    export_data = combine_genetic_data(export_data, 'csv/frontalpole_genes.csv')
    export_all_data(export_data, out_file='csv/PING.csv')
