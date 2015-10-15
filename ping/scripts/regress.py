"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import os

import csv
import matplotlib.pyplot as plt
import numpy as np

from ..ping.apps import PINGSession
from ..ping.apps.regress import find_one_relationship, skip_key, skip_pairing
from ..ping.utils import do_and_plot_regression
from ..research.asymmetry import is_ai_key
from ..research.data import get_derived_data
from ..research.regress import search_all_vs_itself, search_all_pairwise


def search_all_pairwise(all_data):
    """For each pair of variables, look for a significant regression slope."""
    results = []
    keys = list(set(all_data.keys()) - set(('SubjID',)))
    for ki in range(len(keys)):
        key1 = keys[ki]
        if skip_key(key1):
            continue

        for ii in range(ki + 1, len(keys)):
            key2 = keys[ii]
            if skip_key(key2) or skip_pairing(key1, key2):
                continue

            result = find_one_relationship(all_data, key1, key2, rsq_thresh=0.10)
            if result is not None:
                results.append(result)

    # Now, output the sorted result.
    for key, p, r in sorted(results, lambda v1, v2: int(10000 * (abs(v1[2]) - abs(v2[2])))):
        print("Significant at %.2e (r=%.3f): %s" % (p, r, key))



def search_all_vs_one(all_data, key, rsq_thresh=0., covariates=[], plot=False):
    """For each pair of variables, look for a significant regression slope."""

    results = []
    all_keys = list(set(all_data.keys()) - set(('SubjID', key)))
    for all_key in all_keys:
        if not is_ai_key(all_key):
            continue
        try:
            result = find_one_relationship(all_data,
                                           key,
                                           all_key,
                                           covariates=covariates,
                                           rsq_thresh=rsq_thresh,
                                           plot=plot)
        except TypeError as te:
            # Happens when the data type is string or other non-regressable.
            continue
        if result is not None:
            results.append(result)

    # Now, output the sorted result.
    for all_key, p, r in sorted(results, lambda v1, v2: int(10000 * (v1[2] - v2[2]))):
        print("Significant at %.2e (r=%.3f): %s" % (p, r, all_key))



def print_legend():
    legend = {
        'CGH': "Parahippocampal Cingulum",
        'CST': "Corticospinal Tract",
        'IFO': "Inferior-Fronto-Occipital Fasiculus",
        "IFSFC": "Inferior Frontal Superior Frontal Cortex",
        "SCS": "Superior Corticostriate",
        "SLF": "Superior Longitudinal Fasiculus"}

    for key, val in legend.items():
        print('%5s: %s' % (key, val))



def do_regress(*args):
    local = False

    if local:
        all_data = get_derived_data(prefix=SomeData.IMAGING_PREFIX)
        # search_all_pairwise(all_data)

        # search_all_vs_one(all_data, key='MRI_cort_area_ctx_total_LH_PLUS_RH', rsq_thresh=0.0015, plot=False)
        # search_all_vs_one(all_data, key='MRI_cort_area_TOTAL_AI', rsq_thresh=0.005, plot=True)
        # search_all_vs_one(all_data, key='DTI_fiber_vol_CgH_AI', rsq_thresh=0.005)
        search_all_vs_one(all_data, key='DTI_fiber_vol_CgH_AI',
                          covariates=['DTI_fiber_vol_AllFibnoCC_AI'],
                          rsq_thresh=0.005)
        print_legend()

    else:
        try:
            plt.figure()
        except:
            print("Plotting not available.")
            plot = False
        else:
            print("Plotting detected and will be used!")
            plot = True
            plt.close()

        covariates = ['Age_At_IMGExam', 'Gender', 'FDH_23_Handedness_Prtcpnt']#, 'MRI_cort_area_ctx_total_LH_PLUS_RH']
        cache_dir = 'data/regress'

        # search_all_pairwise(plot=plot, cache_dir=cache_dir, covariates=covariates + ['MRI_cort_area_ctx_total_LH_PLUS_RH'])
        # search_all_vs_one(key='MRI_cort_area_ctx_total_LH_PLUS_RH', plot=plot, cache_dir=cache_dir, covariates=covariates)
        search_all_vs_itself(plot=plot, cache_dir=cache_dir, covariates=covariates, abort_if_done=False)
        print_legend()

        if plot:
            plt.show()

        plt.show()


if __name__ == '__main__':
    raise NotImplementedError("Broke this in the refactor.")
