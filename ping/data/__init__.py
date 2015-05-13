"""
"""
import collections
import copy
import csv
import itertools
import os
import warnings

import numpy as np
import pandas
# import statsmodels.formula.api as smf

from ..apps import PINGSession

raw_input = raw_input if 'raw_input' in dir() else input


anatomical_name = {
    'bankssts': 'superior temporal sulcus',
    'caudalanteriorcingulate': 'caudal anterior cingulate',
    'caudalmiddlefrontal': 'caudal middle frontal',
    'frontalpole': 'frontal pole',
    'inferiorparietal': 'inferior parietal',
    'inferiortemporal': 'inferior temporal',
    'isthmuscingulate': 'isthmus cingulate',
    'lateraloccipital': 'lateral occipital',
    'lateralorbitofrontal': 'lateral orbitofrontal',
    'medialorbitofrontal': 'medial orbitofrontal',
    'middletemporal': 'middle temporal',
    'parahippocampal': 'pars hippocampal',
    'parsopercularis': 'pars opercularis',
    'parsorbitalis': 'pars orbitalis',
    'parstriangularis': 'pars triangularis',
    'posteriorcingulate': 'posterior cingulate',
    'rostralanteriorcingulate': 'rostral anterior cingulate',
    'rostralmiddlefrontal': 'rostral middle frontal',
    'superiorfrontal': 'superior frontal',
    'superiorparietal': 'superior parietal',
    'superiortemporal': 'superior temporal',
    'supramarginal': 'supramarginal',
    'temporalpole': 'temporal pole',
    'transversetemporal': 'transverse temporal',


    'ATR': 'anterior thalamic radiations',
    'CST': 'cortico-spinal',
    'CgC': 'cingulum (cingulate)',
    'CgH': 'cingulum (parahippocampal)',
    'CC': 'corpus callosum',
    'Fx': 'fornix',
    'Fxcut': 'fornix (no fimbria)',
    'IFO': 'inferior-fronto-occipital fasciculus',
    'IFSFC': 'inferior frontal superior frontal cortex',
    'ILF': 'inferior longitudinal fasciculus',
    'SCS': 'superior cortico-striate',
    'SIFC': 'striatal inferior frontal cortex',
    'SLF': 'superior longitudinal fasciculus',
    'Unc': 'uncinate fasciculus',
    'fSCS': 'superior cortico-striate (frontal)',
    'pSCS': 'superior cortico-striate (parietal)',
    'pSLF': 'superior longitudinal fasciculus (parietal)',
    'tSLF': 'superior longitudinal fasciculus (temporal)',}

anatomical_order = np.asarray([
    'frontalpole',
    'superiorfrontal',
    'rostralmiddlefrontal',
    'lateralorbitofrontal',
    'parsorbitalis',
    'parstriangularis',
    'parsopercularis',
    'caudalmiddlefrontal',
    'precentral',
    'postcentral',
    'supramarginal',
    'superiorparietal',
    'inferiorparietal',
    'lateraloccipital',
    'inferiortemporal',
    'middletemporal',
    'superiortemporal',
    'transversetemporal',
    'temporalpole',
    #
    'entorhinal',
    'parahippocampal',
    'fusiform',
    'lingual',
    'pericalcarine',
    'cuneus',
    'precuneus',
    'isthmuscingulate',
    'posteriorcingulate',
    'paracentral',
    'caudalanteriorcingulate',
    'rostralanteriorcingulate',
    'medialorbitofrontal',

    # Fiber tracts
    'ATR', #: 'anterior thalamic radiations',
    'CST', #: 'cortico-spinal',
    'CgC', #: 'cingulum (cingulate)',
    'CgH', #: 'cingulum (parahippocampal)',
    'Fx', #: 'fornix',
    'Fxcut', #: 'fornix (no fimbria)',
    'IFO', #: 'inferior-fronto-occipital fasiculus',
    'IFSFC', #: 'inferior frontal superior frontal cortex',
    'ILF', #: 'inferior longitudinal fasiculus',
    'SCS', #: 'superior cortico-striate',
    'SIFC', #: 'striatal inferior frontal cortex',
    'SLF', #: 'superior longitudinal fasiculus',
    'Unc', #: 'uncinate fasiculus',
    'fSCS', #: 'superior cortico-striate (frontal)',
    'pSCS', #: 'superior cortico-striate (parietal)',
    'pSLF', #: 'superior longitudinal fasiculus (parietal)',
    'tSLF', #: 'superior longitudinal fasiculus (temporal)',}


    # Subcortical
    'Accumbens.area',
    'Amygdala',
    'Hippocampus',
    'Caudate',
    'Pallidum',
    'Putamen',
    'Thalamus.Proper'])


def get_anatomical_name(key):
    global anatomical_name
    return  anatomical_name.get(key, key)


def anatomical_sort(keys):
    """Group keys by prefix, then by brain location."""
    # import numpy as np
    # x = np.array([3,5,7,1,9,8,6,6])
    # y = np.array([2,1,5,10,100,6])
    # 
    # index = np.argsort(x)
    # sorted_x = x[index]
    # sorted_index = np.searchsorted(sorted_x, y)
    # 
    # yindex = np.take(index, sorted_index, mode="clip")
    # mask = x[yindex] != y
    # 
    # result = np.ma.array(yindex, mask=mask)
    # print result    
    
    global anatomical_order

    # Parse keys into their anatomical roots
    keys = sorted(keys)
    prefix_of = lambda k: k.split('.')[0]
    key_prefix = [prefix_of(key) for key in keys]
    anatomical_keys = np.asarray([key[(len(prefix_of(key)) + 1):]
                                  for key in keys])
    print(anatomical_keys)

    # Map across the arrays
    index = np.argsort(anatomical_order)
    sorted_ao = anatomical_order[index]
    sorted_index = np.searchsorted(anatomical_order, anatomical_keys)

    keys_index = np.take(index, sorted_index, mode="clip")
    missing_mask = anatomical_order[keys_index] != anatomical_keys

    result = anatomical_order[keys_index]
    result[missing_mask] = anatomical_keys[missing_mask]
    print(result)

    # Everything else that remains is added alphabetically.
    return result
    

def get_lh_key_from_rh_key(key):
    return key.replace('.rh.', '.lh.').replace('.Right.', '.Left.').replace('.R_', '.L_')


def get_rh_key_from_lh_key(key):
    return key.replace('.lh.', '.rh.').replace('.Left.', '.Right.').replace('.L_', '.R_')


def get_nonhemi_key(key):
    return get_rh_key_from_lh_key(key) \
        .replace('.rh.', '.').replace('.Right.', '.').replace('.R_', '_')


def is_nonimaging_key(key):
    return not np.any([key.startswith(p) for p in PINGData.IMAGING_PREFIX])


def get_bilateral_hemi_keys(key):
    """ Given one hemisphere's property name,
    return both."""

    if 'Left' in key or 'Right' in key:
        left_key = key.replace('Right', 'Left')
        right_key = key.replace('Left', 'Right')
    elif '.lh.' in key or '.rh.' in key:
        left_key = key.replace('.rh.', '.lh.')
        right_key = key.replace('.lh.', '.rh.')
    elif '.L_' in key or '.R_' in key:
        left_key = key.replace('.R_', '.L_')
        right_key = key.replace('.L_', '.R_')
    else:
        raise ValueError("Unknown format for key='%s'" % key)

    return left_key, right_key


def col2prop(col_name):
    return col_name.replace('-', '.')


def which_hemi(key):
    rh_markers = ['.Right.', '.rh.', '.R_']
    lh_markers = ['.Left.', '.lh.', '.L_']
    if np.any([m in key for m in rh_markers]):
        return 'rh'
    elif np.any([m in key for m in lh_markers]):
        return 'lh'
    else:
        return None


def filter_data(data_dict, filter_fns, tag=None):
    """Filter functions must take a key and ndarray of values.
    They can return a single boolean (False) to accept/reject the values,
    or an ndarray to select the indices.

    Filters using 'and' logic.
    """
    if filter_fns is None or len(data_dict) == 0:
        return data_dict
    elif not isinstance(filter_fns, dict):
        filter_fns = dict(zip(data_dict.keys(), itertools.cycle([filter_fns])))

    # Limit based on group value and filtering function
    col_len = len(next(iter(data_dict.values())))
    selected_idx = np.ones((col_len,), dtype=bool)
    selected_keys = list(data_dict.keys())
    for filter_key, filter_fn in filter_fns.items():
        if filter_key not in data_dict:
            warnings.warn('%s has a filter, but not found in data; skipping.' % filter_key)
            continue

        # For each key, can have multiple filters
        if not isinstance(filter_fn, collections.Iterable):
            filter_fn = [filter_fn]

        for fn in filter_fn:
            # If we got here, we have a valid key and a valid fn
            filter_vals = data_dict[filter_key]
            filter_output = fn(filter_key, filter_vals[selected_idx])
            if isinstance(filter_output, bool) or isinstance(filter_output, np.bool_):
                if not filter_output:
                    selected_keys.remove(filter_key)
                break  # accepted as whole, skip below logic for efficiency.
            new_idx = np.zeros(selected_idx.shape)
            new_idx[selected_idx] = filter_output
            selected_idx = np.logical_and(selected_idx, new_idx)
            assert  np.any(selected_idx), 'All items were filtered out.'

    # Pull out filtered values and 'tag' the key
    # (this makes documenting results much easier)
    filtered_data = dict()
    for key in selected_keys:
        if key in ['SubjID'] or tag is None:
            tagged_key = key
        else:
            tagged_key = '%s_%s' % (key, tag)

        filtered_data[tagged_key] = np.asarray(data_dict[key])[selected_idx]

    return filtered_data


def merge_by_key(dict1, dict2, merge_key='SubjID', tag=None):

    # Make index mapping from dict1 into dict2
    all_idx = []
    for val in dict1[merge_key]:
        if val in dict2[merge_key]:
            all_idx.append(dict2[merge_key].tolist().index(val))
        else:
            all_idx.append(np.nan)
    all_idx = np.asarray(all_idx)

    # Now reorder dict2 values to dict1
    out_dict = copy.deepcopy(dict1)
    for key, vals in dict2.items():
        # Don't reassign primary key
        if key == merge_key:
            continue

        # Reorder dict2's values into dict1's order
        if np.any(np.isnan(all_idx)):
            reordered_vals = np.asarray([vals[idx] if not np.isnan(idx) else np.nan
                                         for idx in all_idx])
        else:
            reordered_vals = vals[all_idx]

        # Assign, possibly after tagging keys
        if tag is not None:
            key = tag % key
        if key in out_dict:
            pass  # warnings.warn('Overwriting "%s" values from original dict in merge.' % key)
        out_dict[key] = reordered_vals

    return out_dict


class PINGData(object):
    """Dictionary representing a CSV database. Underlying storage is
    a dictionary of parallel arrays, with "SubjID" representing the primary key.
    """
    PING_DATA = None  # shared
    IMAGING_PREFIX = ['MRI_cort_area.ctx', 'MRI_cort_thick.ctx',
                      'MRI_subcort_vol', 'DTI_fiber_vol',
                      'DTI_fiber_FA', 'DTI_fiber_LD', 'DTI_fiber_TD']

    def __init__(self, data=None, scrub_keys=False, scrub_values=True, csv_path=None, username=None, passwd=None, force=False):
        if data is not None:
            pass

        elif PINGData.PING_DATA is not None and not force:
            data = PINGData.PING_DATA

        else:
            csv_path = csv_path or os.path.join('csv', 'PING_raw_data.csv')

            # Download data
            sess = PINGSession(username=username, passwd=passwd)
            if not os.path.exists(csv_path):
                print("Downloading PING data...")
                sess.login()
                sess.download_PING_spreadsheet(out_file=csv_path)
            sess.clean_PING_spreadsheet(out_file=csv_path)

            print("Loading PING data...")
            try:
                data = pandas.read_csv(csv_path, low_memory=False)
            except ValueError as ve:
                # Corrupt spreadsheet. Re-GET
                print("Error loading the PING data: %s" % ve)
                yn = raw_input("The PING spreadsheet is corrupt. Delete and download? (y/N) > ")
                if yn.lower() == 'y':
                    os.remove(csv_path)
                    warnings.warn("How to do recursive calls in OOP... in __init__???")
                    self.__init__(self, scrub_keys=scrub_keys, scrub_values=scrub_values,
                                  csv_path=csv_path, username=username, passwd=passwd,
                                  force=True)
                    return
                raise ve

            # Convert dots to underscores
            print("Converting PING data...")
            new_data = dict()
            for key, val in data.iteritems():
                if scrub_keys:
                    key = key.replace('.', '_')
                if scrub_values:
                    val = val.as_matrix()
                new_data[key] = val
            data = new_data

            # print("Regressing data on confounds...")
            # for key in data.keys():
            #     formula = ('%s ~ FDH_Highest_Education + FDH_3_Household_Income +'
            #                '     DeviceSerialNumber + GAF_africa + GAF_amerind +'
            #                '     GAF_eastAsia + GAF_oceania + GAF_centralAsia') % key.replace('.', '_')
            #     try:
            #         resid = smf.ols(formula, data=data).fit().resid
            #         data[key] = resid
            #     except Exception as e:
            #         print("Failed (%s): %s" % (key, e))
            PINGData.PING_DATA = data

        self.data_dict = copy.deepcopy(data)

    def merge(self, data_dict, merge_key='SubjID', tag=None):
        if isinstance(data_dict, PINGData):
            data_dict = data_dict.data_dict
        self.data_dict = merge_by_key(self.data_dict, data_dict,
                                      merge_key=merge_key, tag=tag)
        return self
        
    def filter(self, filter_fns, tag=None, op='or'):
        if filter_fns is None or not filter_fns:
            return self
        elif op == 'and':
            self.data_dict = filter_data(self.data_dict, filter_fns=filter_fns,
                                         tag=tag)
        elif op == 'or':
            if not isinstance(filter_fns, collections.Iterable):
                filter_fns = [filter_fns]
            
            old_dict = self.data_dict
            new_dict = dict(SubjID=old_dict['SubjID'])
            for fn in filter_fns:
                new_dict.update(filter_data(old_dict, filter_fns=fn, tag=tag))
            self.data_dict = new_dict

        return self

    def export(self, out_file, partial=True):
        keys = list(self.data_dict.keys())
        if partial:
            keys = list(set(keys) - set(self.PING_DATA.keys()))
            keys += ['SubjID']

        dir_path = os.path.dirname(out_file)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        # Now export as csv
        n_subj = len(next(iter(self.data_dict.values())))
        with open(out_file, 'w') as fp:
            w = csv.writer(fp)
            w.writerow(keys)
            for row_idx in range(n_subj):
                row = []
                for key in keys:
                    row.append(self.data_dict[key][row_idx])
                w.writerow(row)

    def get_num_subjects(self):
        if len(self.data_dict) == 0:
            return np.nan
        else:
            return len(next(iter(self.data_dict.values())))

    def get_tbx_data(self):
        return filter_data(self.data_dict,
                           filter_fns=[ lambda k, v: k.startswith('TBX_'),
                                        lambda k, v: v.dtype.name not in ['string', 'object'] ])

    def get_fdh_data(self):
        return filter_data(self.data_dict,
                           filter_fns=[ lambda k, v: k.startswith('FDH'),
                                        lambda k, v: v.dtype.name not in ['string', 'object'] ])

    def get_twohemi_keys(self, filter_fns=None):
        """Given a key prefix, get all keys that have left/right pairs."""

        data_dict = filter_data(self.data_dict, filter_fns)

        # RH only
        data_dict = filter_data(data_dict, lambda k, v: which_hemi(k) == 'rh')

        # No ventricles
        data_dict = filter_data(data_dict, lambda k, v: 'vent' not in k.lower())

        return np.asarray(list(data_dict.keys()))
