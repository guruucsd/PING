"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from ping.access import col2prop, load_PING_data, get_twohemi_keys, get_asymmetry_index
from ping.utils import do_and_plot_regression


def make_groups(data, grouping_prop_names):
    """ Group data by the unique values of `grouping_prop_name`"""

    if not isinstance(grouping_prop_names, list):
        grouping_prop_names = [grouping_prop_names]

    # Group
    grouping_index = group_names = []
    for gpn in grouping_prop_names:
        grouping_data = np.asarray(data[gpn].tolist())
        prop_groups = np.asarray(list(set(np.unique(grouping_data)) -
                                 set(['Not yet established', 'nan'])))
        n_prop_group = len(prop_groups)

        if len(group_names) == 0:
            group_names = prop_groups.tolist()
            grouping_index = [grouping_data == pg for pg in prop_groups]
        else:
            group_names = ['%s_%s' % (g1, pg)
                           for g1 in group_names
                           for pg in prop_groups]
            grouping_index = [np.logical_and(idx1, grouping_data == pg)
                              for idx1 in grouping_index
                              for pg in prop_groups]

    return group_names, grouping_index


def compare_group_asymmetry(data, prop_name, grouping_prop_names):
    """ Groups data according to grouping_prop_names, computes
    asymmetry index for prop_name, and graphs."""

    group_names, grouping_index = make_groups(data, grouping_prop_names)

    age_data = np.asarray(data['Age_At_IMGExam'].tolist())
    fh = plt.figure(figsize=(18, 10))
    n_subplots = 1 + len(group_names)
    n_rows = int(np.round(0.75 * np.sqrt(n_subplots)))
    n_cols = int(np.ceil(n_subplots / float(n_rows)))

    ymeans = []
    ystds = []
    for gi, group_name in enumerate(['all'] + group_names):
        # Index the current group
        if group_name == 'all':
            idx = np.ones((len(data.values()[0]),), dtype=bool)
        else:
            idx = grouping_index[gi - 1]

        # Select data within the group
        cur_ages = age_data[idx]
        prop_asymmetry = get_asymmetry_index(data, prop_name)[idx]

        # Remove bad data
        bad_idx = np.logical_or(np.isnan(cur_ages), np.isnan(prop_asymmetry))
        good_idx = np.logical_not(bad_idx)
        cur_ages = cur_ages[good_idx]
        prop_asymmetry = prop_asymmetry[good_idx]

        # Plot the result
        ax = fh.add_subplot(n_rows, n_cols, gi + 1)
        ymn, ystd = do_and_plot_regression(cur_ages, prop_asymmetry, ax=ax,
                                           xlabel='Age', ylabel='Asymmetry Index (LH - RH)',
                                           title='Group: %s\n%s' % (group_name, prop_name))


def loop_show_asymmetry(prefix, grouping_prop_names=['FDH_23_Handedness_Prtcpnt', 'Gender']):
    """ Loop over all properties to show asymmetry."""
    data = load_PING_data()

    # Process & plot the data.
    for pi, prop_name in enumerate(get_twohemi_keys(prefix, data.keys())):
        compare_group_asymmetry(data, prop_name=prop_name,
                                grouping_prop_names=grouping_prop_names)
        print pi
    plt.show()


if __name__ == '__main__':
    # Grey volume vs. thickness vs. area
    # vol_name = col2prop('MRI_subcort_vol-Right-Cerebral-Cortex')
    # thick_name = col2prop('MRI_cort_thick-ctx-rh-mean')
    # area_name = col2prop('MRI_cort_area-ctx-rh-total')
    # loop_show_asymmetry(prefix=[vol_name, thick_name, area_name])

    # Results: MRI_subcort_vol-Right-Cerebral-Cortex
    # loop_show_asymmetry(prefix='MRI_subcort_vol')

    # Results: 
    # loop_show_asymmetry(prefix='MRI_cort_thick')

    # Results: 
    # loop_show_asymmetry(prefix='MRI_cort_area')

    # Results: A LOT
    # loop_show_asymmetry(prefix='DTI_fiber_vol')

    # Results: Examine on sex
    # loop_show_asymmetry(prefix='DTI_fiber_vol')

    # Doesn't work... sex is missing!
    # loop_show_asymmetry(prefix='DTI_fiber_vol', grouping_prop_names=['Gender', 'FDH_23_Handedness_Prtcpnt'])

    #
    #loop_show_asymmetry(prefix=['MRI_cort_area', 'MRI_cort_thick', 'MRI_subcort_vol', 'DTI_fiber_vol'],
    #                    grouping_prop_names=['Gender', 'FDH_23_Handedness_Prtcpnt'])

    loop_show_asymmetry(prefix='MRI_cort_area_ctx_rh_frontalpole')
