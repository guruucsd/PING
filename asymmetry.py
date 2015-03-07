"""
File for investigating asymmetry from PING data
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas
from scipy.stats import linregress, pearsonr
from sklearn.decomposition import PCA


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


def make_groups(data, grouping_prop_name):

    if not isinstance(grouping_prop_name, list):
        grouping_prop_name = [grouping_prop_name]

    # Group
    grouping_index = group_names = []
    for gpn in grouping_prop_name:
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


def get_asymmetry_index(data, prop_name):
    if 'Left' in prop_name or 'Right' in prop_name:
        left_prop_name = prop_name.replace('Right', 'Left')
        right_prop_name = prop_name.replace('Left', 'Right')
    elif '.lh.' in prop_name or '.rh.' in prop_name:
        left_prop_name = prop_name.replace('.rh.', '.lh.')
        right_prop_name = prop_name.replace('.lh.', '.rh.')
    elif '.L_' in prop_name or '.R_' in prop_name:
        left_prop_name = prop_name.replace('.R_', '.L_')
        right_prop_name = prop_name.replace('.L_', '.R_')
    else:
        raise ValueError("Unknown format for prop_name='%s'" % prop_name)

    # Select data within the group
    LH_data = np.asarray(data[left_prop_name].tolist())
    RH_data = np.asarray(data[right_prop_name].tolist())

    # Compute an asymmetry index
    prop_asymmetry = asymmetry_index(LH_data, RH_data)
    return prop_asymmetry


def compare_group_asymmetry(data, prop_name, grouping_prop_name):
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
        prop_asymmetry = get_asymmetry_index(data, prop_name)[idx]

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


def load_data():
    # Load data
    print("Loading data...")
    script_dir = os.path.abspath(os.path.dirname(__file__))
    csv_path = os.path.join(script_dir, 'PING_raw_data.csv')
    data = pandas.read_csv(csv_path)
    return data


def matching_keys(prefix, all_keys):
    # Massage inputs
    if not isinstance(prefix, list):
        prefix = [prefix]
    prefix = np.asarray(prefix)
    rh_markers = ['Right', '.rh.', '.R_']

    good_keys = []
    for prop_name in all_keys:
        if (not np.any(np.asarray([r in prop_name for r in rh_markers])) or
                not np.any(np.asarray([p in prop_name for p in prefix]))):
            continue
        if 'vent' in prop_name.lower() and 'ventral' not in prop_name.lower():
            print("Skipping %s" % prop_name)
            continue
        good_keys.append(prop_name)
    return np.asarray(good_keys)


def loop_show_asymmetry(prefix, grouping_prop_name='FDH_23_Handedness_Prtcpnt'):

    data = load_data()

    # Process & plot the data.
    for pi, prop_name in enumerate(matching_keys(prefix, data.keys())):
        compare_group_asymmetry(data, prop_name=prop_name,
                                grouping_prop_name=grouping_prop_name)
        pi += 1
        print pi
    plt.show()


def show_asymmetry_by_component(prefix, grouping_prop_name='FDH_23_Handedness_Prtcpnt'):
    """ Grab all data matching prefix, split into groups by grouping_prop_name,
    then perform PCA and show the PCs."""

    data = load_data()
    good_keys = matching_keys(prefix, data.keys())

    # Narrow to selected keys
    data_subset = []
    for key in good_keys:
        data_subset.append(get_asymmetry_index(data, key))
    data_subset = np.asarray(data_subset)

    # Eliminate subjects with zero asymmetry in anything (likely artifacts)
    noasymm_idx = np.abs(data_subset).sum(axis=0) != 0.
    print("Removing %d subjects with ZERO asymmetry "
          "(over %d measures!)" % (int(np.logical_not(noasymm_idx).sum()),
                                   len(good_keys)))
    data_subset = data_subset[:, noasymm_idx]

    # Narrow to selected keys with mean(asymmetry_index) > 0.05
    mean_asymmetry = data_subset.mean(axis=1)
    asymmetric_idx = np.abs(mean_asymmetry) >= 0.05
    data_subset = data_subset[asymmetric_idx]
    good_keys = good_keys[asymmetric_idx]
    mean_asymmetry = mean_asymmetry[asymmetric_idx]

    # Report a summary of the asymmetries
    n_measures = len(good_keys)
    all_thresh = [1.0, 0.25, 0.10, 0.05, 0.01]
    for ti, thresh in enumerate(all_thresh[1:]):
        abs_mean = np.abs(mean_asymmetry)
        idx_above = np.logical_and(thresh <= abs_mean, abs_mean < all_thresh[ti])
        n_above = int(idx_above.astype(int).sum())
        print("%d (of %d) asymmetries > %.2f" % (n_above, n_measures, thresh))
        for hemi_sign, hemi_name in zip([1, -1], ['lh', 'rh']):
            hemi_idx = np.logical_and(np.sign(mean_asymmetry) == hemi_sign, idx_above)
            n_hemi = np.abs(hemi_idx.sum())
            print("\t%d (of %d) are %s" % (n_hemi, n_above, hemi_name))
            for idx in np.nonzero(hemi_idx.astype(int))[0]:
                print("\t\t%s (%.2f)" % (good_keys[idx], mean_asymmetry[idx]))

    # Compute PCA and dump the results
    pca = PCA(whiten=False)
    pca.fit(data_subset.T)
    selected_idx = pca.explained_variance_ratio_ >= 0.05
    for pi, pc in enumerate(pca.components_[selected_idx]):
        print "%2d: (%.2f):" % (pi, pca.explained_variance_ratio_[pi])
        for hemi_sign, hemi_name in zip([1, -1], ['lh', 'rh']):
            for key, coeff, mag in zip(good_keys, pc, mean_asymmetry):
                if np.sign(mag) == hemi_sign:
                    print "\t%s %s%.4f %s (%.2f)" % (hemi_name, ' ' if coeff >= 0 else '', coeff, key, 1000*mag*coeff)

    # Now, take this factor and compare it to behavioral data
    # (or, should I throw the behavior in the PCA?)
    tbx_props = filter(lambda k: 'TBX_' in k, data.keys())
    tbx_data = np.asarray(data[tbx_props].to_records())[noasymm_idx]
    for pi, pc in enumerate(pca.components_[selected_idx]):
        pc_projection = np.dot(pc, data_subset)
        for ti, tbx_prop in enumerate(tbx_props):
            row_data = np.asarray([r[ti] for r in tbx_data])
            if 'string' in row_data.dtype.name:
                continue
            good_subj = np.logical_not(np.isnan(row_data))
            rval, pval = pearsonr(pc_projection[good_subj], row_data[good_subj])
            if pval < 0.05:
                print('%d: %s %.4f, %.4f' % (pi, tbx_prop, rval, pval))


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
    # loop_show_asymmetry(prefix='MRI_cort_thick')

    # Results: 
    # loop_show_asymmetry(prefix='MRI_cort_area')

    # Results: A LOT
    # loop_show_asymmetry(prefix='DTI_fiber_vol')

    # Results: Examine on sex
    # loop_show_asymmetry(prefix='DTI_fiber_vol')

    # Doesn't work... sex is missing!
    #loop_show_asymmetry(prefix='DTI_fiber_vol', grouping_prop_name=['Gender', 'FDH_23_Handedness_Prtcpnt'])

    #
    show_asymmetry_by_component(prefix=['MRI_subcort_vol', 'MRI_cort_thick', 'MRI_cort_area', 'DTI_fiber_vol'])
