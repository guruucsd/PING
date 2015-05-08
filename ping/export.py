"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import copy
import csv
import os

import matplotlib.pyplot as plt
import numpy as np


def merge_by_key(dict1, dict2, merge_key='SubjID'):

    # Make index mapping from dict1 into dict2
    all_idx = []
    for val in dict1[merge_key]:
        if val in dict2[merge_key]:
            all_idx.append(dict2[merge_key].tolist().index(val))
        else:
            all_idx.append(np.nan)

    # Now reorder dict2 values to dict1
    out_dict = copy.copy(dict1)
    for key, vals in dict2.items():
        if key == merge_key:
            continue
        reordered_vals = [vals[idx] if not np.isnan(idx) else np.nan
                          for idx in all_idx]
        out_dict[key] = reordered_vals

    return out_dict


def update_data_dictionary(data_dict, update_dict=None, update_csv=None,
                           data_type=None):
    assert int(update_dict is None) + int(update_csv is None) == 1, "Dict or file"

    if update_dict is not None:
        if 'SubjID' in data_dict and 'SubjID' in update_dict:
            out_dict = merge_by_key(data_dict, update_dict, merge_key='SubjID')
        elif len(data_dict.values()[0]) == len(update_dict.values()[0]):
            # Same length, assume they match.
            out_dict = copy.copy(data_dict)
            out_dict.update(update_dict)
        else:
            raise ValueError('Cannot combine dicts; no SubjID field and data are different sizes.')

    else:
        # Below here, using csv_file
        if data_type is None:
            # Merge the two together.
            with open(gene_file, 'rb') as fp:
                csv_data = np.recfromcsv(fp)
            out_dict = copy.copy(data_dict)
            out_dict.update(csv_data)
        elif data_type == 'genetic':
            out_dict = copy.copy(data_dict)
            out_dict = combine_genetic_data(out_dict, csv_file)
        else:
            raise NotImplementedError(data_type)

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
