"""
"""
import copy
import os

import numpy as np
import pandas

from ping.access import (get_lh_prop_name, get_nonhemi_prop_name,
                         load_PING_data, get_twohemi_keys)
from ping.asymmetry import (get_asymmetry_index, get_ai_prop_name,
                            is_ai_prop_name)
from ping.export import export_all_data, merge_by_key
from ping.multivariate import AsymmetryPCA

EXPORTED_PING_SPREADSHEET = 'csv/PING_userdefined.csv'


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


def compute_all_asymmetries(prefix):
    """ Loop over all properties to show asymmetry."""
    data = load_PING_data(scrub_fields=False)

    export_data = dict((('SubjID', data['SubjID'],),))

    # Process & plot the data.
    for prop_name in get_twohemi_keys(data.keys(), prefix=prefix):
        dest_prop_name = get_ai_prop_name(prop_name)
        export_data[dest_prop_name] = get_asymmetry_index(data, prop_name,
                                                          mask_nan=False)

    # Add one property for total asymmetry
    n_subj = len(data['SubjID'])
    for p in prefix:
        total_asymmetry = np.zeros((n_subj,))
        good_keys = filter(lambda k: k.startswith(p) and is_ai_prop_name(k),
                           export_data.keys())
        for key in good_keys:
            values = export_data[key].copy()
            values[np.isnan(values)] = 0.
            total_asymmetry += export_data[key]**2
        total_ai_key = get_ai_prop_name(p + '_TOTAL')
        export_data[total_ai_key] = np.sqrt(total_asymmetry)
        assert len(export_data[total_ai_key]) == n_subj

    return export_data


def compute_component_loadings(prefix):
    all_data = compute_all_asymmetries(prefix)

    pca = AsymmetryPCA(whiten=True)
    pca.fit(all_data)

    pca.report_asymmetry_loadings()
    pca.report_behavior_correlations()
    pca.report_background_correlations()

    pc_projections = pca.get_projections()
    keys = ['PC%d' % pi for pi in range(pc_projections.shape[0])]

    pc_dict = dict(zip(keys, pc_projections))
    pc_dict['SubjID'] = pca.subj_ids

    return pc_dict


def combine_genetic_data(export_data, gene_file):
    # Merge the two together.
    with open(gene_file, 'rb') as fp:
        gene_data = np.recfromcsv(fp)

    all_idx = []
    for subj_id in export_data['SubjID']:
        x_subj_id = '"%s"' % subj_id
        if x_subj_id in gene_data['subjid']:
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
                new_data[key.replace('.', '_')] = data[key].as_matrix()
            data = new_data

    if data is None:
        print "Computing derived data..."
        data = compute_all_asymmetries(prefix=prefix)
        data = merge_by_key(data, compute_all_totals(prefix=prefix))

        # Add principle components as a whole,
        #   then for each prefix seperately.
        try:
            pc_dict = compute_component_loadings(prefix=prefix)
        except Exception as e:
            print "Skipping PCA: %s" % e
        else:
            recoded_keys = ['ALL_%s' % k if k != 'SubjID' else k
                            for k in pc_dict.keys()]
            data = merge_by_key(data, dict(zip(recoded_keys, pc_dict.values())))

            for p in prefix:
                pc_dict = compute_component_loadings(prefix=p)
                recoded_keys = ['%s_%s' % (p, k) if k != 'SubjID' else k
                                for k in pc_dict.keys()]
                data = merge_by_key(data, dict(zip(recoded_keys, pc_dict.values())))

        data = combine_genetic_data(data, 'csv/frontalpole_genes.csv')

    return data


def get_all_data(prefix, force=False):
    ping_data = copy.deepcopy(load_PING_data())
    all_data = get_derived_data(prefix=prefix, force=force)

    return merge_by_key(ping_data, all_data)
