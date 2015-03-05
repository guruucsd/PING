"""
File for investigating asymmetry from PING data
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas
from scipy.stats import linregress


def asymmetry_index(left, right):
    aidx = (left - right) / ((left + right) / 2.)
    aidx[np.isnan(aidx)] = 0
    return aidx


def do_and_plot_regression(X, Y, xlabel='', ylabel='', title='', ax=None):
    assert not np.any(np.isnan(X))
    assert not np.any(np.isnan(Y))

    m, b, rval,  pval, stderr = linregress(X, Y)
    if not ax:
        fh = plt.figure()
        ax = ax
    xlim = [X.size and X.min() or 0, X.size and X.max() or 1]
    xvals = np.linspace(xlim[0], xlim[1], 25)
    ax.plot(xvals, m * xvals + b, 'r', linewidth=3.)
    ax.set_title('%s\n(r=%.3f, p=%.3e, n=%d)' % (title, rval, pval, X.size))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.hold('on')
    ax.plot(xvals, xvals * 0., 'k--', linewidth=5)
    ax.scatter(X, Y)
    ax.set_ylim([-1, 1])


def compare_group_asymmetry(data, prop_name, grouping_prop_name):
    if 'Left' in prop_name or 'Right' in prop_name:
        left_prop_name = prop_name.replace('Right', 'Left')
        right_prop_name = prop_name.replace('Left', 'Right')
        ylabel = left_prop_name.replace('Left', '')
    elif '.lh.' in prop_name or '.rh.' in prop_name:
        left_prop_name = prop_name.replace('.rh.', '.lh.')
        right_prop_name = prop_name.replace('.lh.', '.rh.')
        ylabel = left_prop_name.replace('.lh.', '.')
    elif '.L_' in prop_name or '.R_' in prop_name:
        left_prop_name = prop_name.replace('.R_', '.L_')
        right_prop_name = prop_name.replace('.L_', '.R_')
        ylabel = left_prop_name.replace('.L_', '.')
    else:
        raise ValueError("Unknown format for prop_name='%s'" % prop_name)

    grouping_data = np.asarray(data[grouping_prop_name].tolist())
    groups = np.asarray(list(set(np.unique(grouping_data)) - set(['Not yet established', 'nan'])))
    n_subplots = 1 + len(groups)

    age_data = np.asarray(data['Age_At_IMGExam'].tolist())
    fh = plt.figure(figsize=(18, 10))
    for gi, group_name in enumerate(['all'] + groups.tolist()):
        # Index the current group
        if group_name == 'all':
            idx = grouping_data != '__magic_string__'
        else:
            idx = grouping_data == group_name

        # Select data within the group
        cur_ages = age_data[idx]
        LH_data = np.asarray(data[left_prop_name].tolist())[idx]
        RH_data = np.asarray(data[right_prop_name].tolist())[idx]

        # Remove bad data
        bad_idx = np.logical_or(np.isnan(cur_ages), np.isnan(LH_data), np.isnan(RH_data))
        good_idx = np.logical_not(bad_idx)
        cur_ages = cur_ages[good_idx]
        LH_data = LH_data[good_idx]
        RH_data = RH_data[good_idx]

        # Compute an asymmetry index
        prop_asymmetry = asymmetry_index(LH_data, RH_data)

        # Plot the result
        ax = fh.add_subplot(1, n_subplots, gi + 1)
        do_and_plot_regression(cur_ages, prop_asymmetry, ax=ax,
                               xlabel='Age', ylabel='Asymmetry Index (LH - RH)',
                               title='Group: %s\n%s' % (group_name, prop_name))


def loop_show_asymmetry(prefix, grouping_prop_name='FDH_23_Handedness_Prtcpnt'):
    # Massage inputs
    if not isinstance(prefix, list):
        prefix = [prefix]
    prefix = np.asarray(prefix)
    rh_markers = ['Right', '.rh.', '.R_']

    # Load data
    print("Loading data...")
    script_dir = os.path.abspath(os.path.dirname(__file__))
    csv_path = os.path.join(script_dir, 'PING_raw_data.csv')
    data = pandas.read_csv(csv_path)

    # Process & plot the data.
    pi = 0
    for prop_name in data.keys():
        if (not np.any(np.asarray([r in prop_name for r in rh_markers])) or
                not np.any(np.asarray([p in prop_name for p in prefix]))):
            continue
        compare_group_asymmetry(data, prop_name=prop_name,
                                grouping_prop_name=grouping_prop_name)
        pi += 1
        print pi
    plt.show()


def col2prop(col_name):
    return col_name.replace('-', '.')


if __name__ == '__main__':
    # Grey volume vs. thickness vs. area
    vol_name = col2prop('MRI_subcort_vol-Right-Cerebral-Cortex')
    thick_name = col2prop('MRI_cort_thick-ctx-rh-mean')
    area_name = col2prop('MRI_cort_area-ctx-rh-total')
    loop_show_asymmetry(prefix=[vol_name, thick_name, area_name])

    # Results: MRI_subcort_vol-Right-Cerebral-Cortex
    # loop_show_asymmetry(prefix='MRI_subcort_vol')

    # Results: 
    loop_show_asymmetry(prefix='MRI_cort_thick')

    # Results: 
    loop_show_asymmetry(prefix='MRI_cort_area')

    # Results: A LOT
    # loop_show_asymmetry(prefix='DTI_fiber_vol')
