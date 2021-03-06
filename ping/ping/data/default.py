"""
"""
import collections
import copy
import csv
import os

import numpy as np
import pandas
# import statsmodels.formula.api as smf

from .base import filter_data, merge_by_key
from ..apps import PINGSession

# for py2/3 compatibility
raw_input = raw_input if 'raw_input' in dir() else input


class PINGData(object):
    """Dictionary representing a CSV database. Underlying storage is
    a dictionary of parallel arrays, with "SubjID" representing the primary key.
    """
    PING_DATA = None  # shared
    IMAGING_PREFIX = ['MRI_cort_area.ctx', 'MRI_cort_thick.ctx', 'MRI_cort_vol.ctx',
                      'MRI_subcort_vol', 'DTI_fiber_vol',
                      'DTI_fiber_FA', 'DTI_fiber_LD', 'DTI_fiber_TD']

    def __init__(self, data=None, scrub_keys=False, scrub_values=True,
                 csv_path=None, username=None, passwd=None, force=False,
                 data_dir='data'):

        if data is not None:
            # User specified data directly.
            pass

        elif self.PING_DATA is not None and not force:
            # Data has already been loaded.
            data = self.PING_DATA

        else:
            # Get the PING raw data (download if necessary)
            csv_path = csv_path or os.path.join(data_dir, 'PING_raw_data.csv')

            # Download data
            sess = PINGSession(username=username, passwd=passwd, data_dir=data_dir)
            if not os.path.exists(csv_path):
                print("Downloading PING data...")
                sess.login()
                sess.download_PING_spreadsheet(out_file=csv_path)
            sess.clean_PING_spreadsheet(out_file=csv_path)

            print("Loading PING data...")
            try:
                data = pandas.read_csv(csv_path, low_memory=False)
            except ValueError as ve:
                # Corrupt spreadsheet. Delete and re-download
                print("Error loading the PING data: %s" % ve)
                yn = raw_input("The PING spreadsheet is corrupt. Delete and download? (y/N) > ")
                if yn.lower() == 'y':
                    os.remove(csv_path)
                    # warnings.warn("How to do recursive calls in OOP... in __init__???")
                    self.__init__(self, scrub_keys=scrub_keys, scrub_values=scrub_values,
                                  csv_path=csv_path, username=username, passwd=passwd,
                                  force=True)
                    return
                raise ve

            # Convert dots to underscores
            print("Converting PING data...")
            new_data = dict()
            for key in data:
                val = data[key]
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
            self.PING_DATA = data

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

    def purge_empty_subjects(self):
        # Find # of nan numeric measures
        any_good = np.zeros((self.get_num_subjects(),))
        for key in self.data_dict.keys():
            try:
                any_good += ~np.isnan(self.data_dict[key])
            except TypeError:
                continue
        # Remove subjects with no non-nan numeric measures.
        for key in self.data_dict.keys():
            self.data_dict[key] = self.data_dict[key][any_good > 0]
        assert len(next(iter(self.data_dict.values()))) > 0

    def get_num_subjects(self):
        if len(self.data_dict) == 0:
            return np.nan
        else:
            return len(next(iter(self.data_dict.values())))

    def get_tbx_data(self):
        return filter_data(self.data_dict,
                            filter_fns=[lambda k, v: k.startswith('TBX_'),
                                        lambda k, v: v.dtype.name not in ['string', 'object']])

    def get_fdh_data(self):
        return filter_data(self.data_dict,
                           filter_fns=[lambda k, v: k.startswith('FDH'),
                                       lambda k, v: v.dtype.name not in ['string', 'object']])

    def get_twohemi_keys(self, filter_fns=None):
        """Given a key prefix, get all keys that have left/right pairs."""
        data_dict = filter_data(self.data_dict, filter_fns)

        # RH only
        data_dict = filter_data(data_dict, lambda k, v: self.which_hemi(k) == 'rh')

        # No ventricles
        data_dict = filter_data(data_dict, lambda k, v: 'vent' not in k.lower())

        return np.asarray(list(data_dict.keys()))

    anatomical_name = {
        'bankssts': 'superior temporal sulcus',
        'caudalanteriorcingulate': 'caudal anterior cingulate',
        'caudalmiddlefrontal': 'caudal middle frontal',
        'cuneus': 'cuneus',
        'entorhinal': 'entorhinal',
        'fusiform': 'fusiform',
        'frontalpole': 'frontal pole',
        'inferiorparietal': 'inferior parietal',
        'inferiortemporal': 'inferior temporal',
        'isthmuscingulate': 'isthmus cingulate',
        'lateraloccipital': 'lateral occipital',
        'lateralorbitofrontal': 'lateral orbitofrontal',
        'lingual': 'lingual',
        'medialorbitofrontal': 'medial orbitofrontal',
        'middletemporal': 'middle temporal',
        'paracentral': 'paracentral',
        'parahippocampal': 'parahippocampal',
        'parsopercularis': 'pars opercularis',
        'parsorbitalis': 'pars orbitalis',
        'parstriangularis': 'pars triangularis',
        'pericalcarine': 'pericalcarine',
        'postcentral': 'postcentral',
        'posteriorcingulate': 'posterior cingulate',
        'precentral': 'precentral',
        'precuneus': 'precuneus',
        'rostralanteriorcingulate': 'rostral anterior cingulate',
        'rostralmiddlefrontal': 'rostral middle frontal',
        'superiorfrontal': 'superior frontal',
        'superiorparietal': 'superior parietal',
        'superiortemporal': 'superior temporal',
        'supramarginal': 'supramarginal',
        'temporalpole': 'temporal pole',
        'transversetemporal': 'transverse temporal',

        'IFSFC': 'inferior frontal superior frontal cortex',
        'SIFC': 'striatal inferior frontal cortex',
        'Unc': 'uncinate fasciculus',
        'ATR': 'anterior thalamic radiations',
        'fSCS': 'superior cortico-striate (frontal)',
        'IFO': 'inferior-fronto-occipital fasciculus',
        'CgC': 'cingulum (cingulate)',
        'ILF': 'inferior longitudinal fasciculus',
        'SCS': 'superior cortico-striate',
        'pSCS': 'superior cortico-striate (parietal)',
        'pSLF': 'superior longitudinal fasciculus (parietal)',
        'tSLF': 'superior longitudinal fasciculus (temporal)',
        'SLF': 'superior longitudinal fasciculus',
        'CgH': 'cingulum (parahippocampal)',
        'CST': 'cortico-spinal',
        'CC': 'corpus callosum',
        'Fx': 'fornix',
        'Fxcut': 'fornix (no fimbria)'}

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
        'IFSFC', #: 'inferior frontal superior frontal cortex',
        'SIFC', #: 'striatal inferior frontal cortex',
        'Unc', #: 'uncinate fasiculus',
        'ATR', #: 'anterior thalamic radiations',
        'fSCS', #: 'superior cortico-striate (frontal)',
        'IFO', #: 'inferior-fronto-occipital fasiculus',
        'CgC', #: 'cingulum (cingulate)',
        'ILF', #: 'inferior longitudinal fasiculus',
        'SCS', #: 'superior cortico-striate',
        'pSCS', #: 'superior cortico-striate (parietal)',
        'pSLF', #: 'superior longitudinal fasiculus (parietal)',
        'tSLF', #: 'superior longitudinal fasiculus (temporal)',}
        'SLF', #: 'superior longitudinal fasiculus',
        'CgH', #: 'cingulum (parahippocampal)',
        'CST', #: 'cortico-spinal',
        'CC', #: 'cortico-spinal',
        'Fx', #: 'fornix',
        'Fxcut', #: 'fornix (no fimbria)',

        # Subcortical
        'Accumbens.area',
        'Amygdala',
        'Hippocampus',
        'Caudate',
        'Pallidum',
        'Putamen',
        'Thalamus.Proper'])

    @classmethod
    def get_anatomical_name(klass, key):
        """Returns a proper anatomical name, or the key if not found."""
        try:
            normd_key = klass.norm_key(key)
            anat_key = klass.anatomical_name.get(normd_key, normd_key)
            if '_' in str(anat_key):
                # no change
                pass  # import pdb; pdb.set_trace()
            return anat_key
        except AttributeError as ae:
            import pdb; pdb.set_trace()
            return None

    @classmethod
    def norm_key(klass, key):
        return klass.norm_keys([key])[0]

    @classmethod
    def get_prefix(klass, key):
        return klass.get_prefixes([key])[0]

    @classmethod
    def prefix2text(klass, prefix):
        d = {
            'MRI_cort_thick.ctx.': 'Cortical thickness (mm)',
            'MRI_cort_area.ctx.': 'Cortical surface area (mm^2)',
            'MRI_subcort_vol.': 'Subcortical volume (mm^3)',
            'DTI_fiber_vol.': 'Fiber tract volume (via DTI) (mm^3)',
            'DTI_fiber_FA.': 'Fiber tract fractional anisotropy (FA)'}

        if prefix in d:
            return d[prefix]
        imperfect_matches = filter(lambda k: k.startswith(prefix), d.keys())
        if imperfect_matches:
            return d[imperfect_matches[0]]
        return prefix

    @classmethod
    def get_measure_key(klass, common_key, measure_keys):
        """First hit in an array (none otherwise)"""
        hits = [k for k in measure_keys if common_key in k]
        if len(hits) == 0:
            return None
        else:
            return hits[0]

    @classmethod
    def get_prefixes(klass, keys):
        normd_keys = sorted([klass.get_nonhemi_key(k) for k in keys])

        # Select every prefix match with the key
        # Fair warning: I add a dummy character to avoid
        #   adding one on the prefix in the next step...
        #   which would be tough with the empty string logic!
        # Get the substring without the prefix (or the key itself
        #   if no prefix was found)
        key_prefixes = [[''] + [p + "." for p in klass.IMAGING_PREFIX
                                if k.startswith(p)]
                        for k in normd_keys]
        return key_prefixes

    @classmethod
    def norm_keys(klass, keys):
        """remove prefix, get nonhemi version"""

        normd_keys = sorted([klass.get_nonhemi_key(k) for k in keys])
        key_prefixes = klass.get_prefixes(normd_keys)
        anatomical_keys = np.asarray([k[len(p[-1]):]
                                      for k, p in zip(normd_keys, key_prefixes)])

        assert len(keys) == len(anatomical_keys), "Return as many keys as requested."
        assert np.all([len(k) > 0 for k in anatomical_keys]), "No keys should be blank."
        return anatomical_keys

    @classmethod
    def anatomical_sort(klass, keys, regroup_results=True):
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

        # Normalize the form of the key
        keys = np.asarray(keys)
        anatomical_keys = klass.norm_keys(keys)

        # Find map from the structures in anatomical_order into the 'keys' list.
        index = np.argsort(anatomical_keys)
        sorted_ak = anatomical_keys[index]
        sorted_index = np.searchsorted(sorted_ak, klass.anatomical_order)

        # Now use that map to reorder the keys themselves.
        keys_index = np.take(index, sorted_index, mode="clip")
        good_mask = anatomical_keys[keys_index] == klass.anatomical_order

        # Grab the reordered keys, append any keys not found.
        result = keys[keys_index[good_mask]]

        # Find which keys remain, then append.
        all_keys_index = set(np.arange(len(keys)))
        found_keys_index = set(np.unique(keys_index[good_mask]))
        missing_keys_idx = np.asarray(list(all_keys_index - found_keys_index), dtype=int)
        result = np.concatenate([result, keys[missing_keys_idx]])

        if regroup_results:
            # Regroup keys by prefix; useful when there are multiple prefixes.
            regrouped_results = []
            for p in klass.IMAGING_PREFIX:
                regrouped_results += [k for k in result if k.startswith(p)]
            regrouped_results += [k for k in result if k not in regrouped_results]
            result = regrouped_results

        assert len(result) == len(keys)

        # Everything else that remains is added alphabetically.
        return result

    @classmethod
    def get_lh_key_from_rh_key(klass, key):
        return key.replace('.rh.', '.lh.').replace('.Right.', '.Left.').replace('.R_', '.L_')

    @classmethod
    def get_rh_key_from_lh_key(klass, key):
        return key.replace('.lh.', '.rh.').replace('.Left.', '.Right.').replace('.L_', '.R_')

    @classmethod
    def get_nonhemi_key(klass, key):
        new_key = klass.get_rh_key_from_lh_key(key) \
                      .replace('.rh.', '.').replace('.Right.', '.').replace('.R_', '_') \
                      .replace('_AI', '').replace('_LH_PLUS_RH', '').replace('_TOTAL', '')  # hacks... for now
        return new_key

    @classmethod
    def is_nonimaging_key(klass, key):
        return not np.any([key.startswith(p) for p in klass.IMAGING_PREFIX])

    @classmethod
    def get_bilateral_hemi_keys(klass, key):
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

    # @classmethod
    # def col2prop(klass, col_name):
    #     return col_name.replace('-', '.')

    @classmethod
    def which_hemi(klass, key):
        rh_markers = ['.Right.', '.rh.', '.R_']
        lh_markers = ['.Left.', '.lh.', '.L_']
        if np.any([m in key for m in rh_markers]):
            return 'rh'
        elif np.any([m in key for m in lh_markers]):
            return 'lh'
        else:
            return None

    @classmethod
    def get_roi_key_from_name(klass, anatomical_name):
        for key, val in klass.anatomical_name.items():
            if val.startswith(anatomical_name):  # allow partial matches...
                return key
        return None
