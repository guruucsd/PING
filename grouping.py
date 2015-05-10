"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy

from ping.access import col2prop, load_PING_data, get_twohemi_keys
from ping.asymmetry import get_asymmetry_index
from ping.utils import do_and_plot_regression


def make_groups(data, grouping_prop_names):
    """ Group data by the unique values of `grouping_prop_name`"""

    if not isinstance(grouping_prop_names, list):
        grouping_prop_names = [grouping_prop_names]

    # Group
    grouping_index = group_names = []
    for gpn in grouping_prop_names:
        grouping_data = np.asarray(data[gpn].tolist())
        prop_groups = (set(np.unique(grouping_data)) -
                       set(['Not yet established', 'nan']))
        prop_groups = np.asarray(list(prop_groups))
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


def compare_group_asymmetry(data, prop_name, grouping_prop_names, plots):
    """ Groups data according to grouping_prop_names, computes
    asymmetry index for prop_name, and graphs."""

    group_names, grouping_index = make_groups(data, grouping_prop_names)

    age_data = np.asarray(data['Age_At_IMGExam'].tolist())
    prop_ai = get_asymmetry_index(data, prop_name)

    n_subplots = len(group_names)
    n_rows = int(np.round(np.sqrt(n_subplots)))
    n_cols = int(np.ceil(n_subplots / float(n_rows)))

    if 'regressions' in plots:
        fh1 = plt.figure(figsize=(18, 10))
        fh1.suptitle('%s regression' % prop_name)
    if 'distributions' in plots:
        fh2 = plt.figure(figsize=(18, 10))
        fh2.suptitle('%s distributions' % prop_name)
        n_bins = 15
        bins = np.linspace(prop_ai[~np.isnan(prop_ai)].min(),
                           prop_ai[~np.isnan(prop_ai)].max(),
                           n_bins + 2)[1:-1]

    group_samples = []
    regressions = []
    for gi, group_name in enumerate(group_names):

        # Index the current group
        if group_name == 'all':
            idx = np.ones((len(data.values()[0]),), dtype=bool)
        else:
            idx = grouping_index[gi]

        # Select data within the group
        cur_ages = age_data[idx]
        group_ai = prop_ai[idx]

        # Remove bad data
        bad_idx = np.logical_or(np.isnan(cur_ages), np.isnan(group_ai))
        good_idx = np.logical_not(bad_idx)
        cur_ages = cur_ages[good_idx]
        group_ai = group_ai[good_idx]

        group_samples.append(group_ai)
        regressions.append(scipy.stats.linregress(cur_ages, group_ai))

        # Plot the regression result
        params = dict(xlabel='Age', ylabel='Asymmetry Index (LH - RH)',
                      title='Group: %s (n=%d)' % (group_name, len(group_ai)))

        if 'regressions' in plots:
            ax1 = fh1.add_subplot(n_rows, n_cols, gi + 1)
            do_and_plot_regression(cur_ages, group_ai, ax=ax1, **params)

        # Plot the distribution result
        if 'distributions' in plots:
            ax2 = fh2.add_subplot(n_rows, n_cols, gi + 1)
            x, _, h = ax2.hist(group_ai, normed=True, bins=bins)
            for item in h:
                item.set_height(item.get_height() / sum(x))

            ax2.set_title(params['title'])
            ax2.set_xlabel(params['ylabel'])
    regressions = np.asarray(regressions)

    # stats[:, n]: 0:mean_stat: 1:mean_pval; 2:std_stat, 3:std_pval
    # shape: n_compares x 4
    stats = np.asarray([(scipy.stats.ttest_ind(gsamps1, gsamps2) +
                         scipy.stats.levene(gsamps1, gsamps2))
                        for gi, gsamps1 in enumerate(group_samples)
                        for gsamps2 in group_samples[gi + 1:]])

    # Test whether variances differ
    if 'stats' in plots:
        dist_mat = scipy.spatial.distance.squareform(stats[:, 2])
        dist_mat = np.ma.masked_where(np.triu(dist_mat == np.nan), dist_mat)
        sig_mat = stats[:, 3] <= (0.05 / stats.shape[0])

        fh3 = plt.figure()
        fh3.suptitle(str(['%.2e' % s for s in stats[:, 3]]))

        ax1 = fh3.add_subplot(1, 2, 1)
        plot_symmetric_matrix_as_triangle(dist_mat, ax=ax1, lbls=group_names)

        ax2 = fh3.add_subplot(1, 2, 2, axisbg=fh3.get_facecolor())
        plot_symmetric_matrix_as_triangle(sig_mat, ax=ax2, lbls=group_names)

    if 'distributions' in plots:
        equalize_xlims(fh2)
        equalize_ylims(fh2)

    return group_names, stats, regressions


def plot_symmetric_matrix_as_triangle(mat, ax=None, lbls=None, vmin=0, vmax=1):
    if ax is None:
        ax = plt.figure().gca()

    ax.set_axis_bgcolor(ax.get_figure().get_facecolor())
    ax.imshow(mat, vmin=vmin, vmax=vmax, interpolation='none')
    ax.set_frame_on(False)

    if lbls:
        sz = mat.shape[0]
        ax.set_xticks(range(sz - 1))
        ax.set_xticklabels(lbls[:sz - 1])
        ax.set_yticks(range(1, sz))
        ax.set_yticklabels(lbls[1:])


def equalize_xlims(fh):
    xlims = np.asarray([ax.get_xlim() for ax in fh.get_axes()])
    xlim = [xlims[:, 0].min(), xlims[:, 1].max()]

    for ax in fh.get_axes():
        ax.set_xlim(xlim)


def equalize_ylims(fh):
    ylims = np.asarray([ax.get_ylim() for ax in fh.get_axes()])
    ylim = [ylims[:, 0].min(), ylims[:, 1].max()]

    for ax in fh.get_axes():
        ax.set_ylim(ylim)


def dump_regressions_csv(regressions, group_names, measure_names):
    # Dump a tsv of group rvals, pvals, and coeff of variation
    n_measures = len(regressions)
    stat_names = ['rval', 'pval', 'coeff_of_var']
    for mi, measure_name in enumerate(sorted(measure_names)):
        if mi == 0:
            header_vals = ['%s_%s' % (group_name, stat_name)
                           for si, stat_name in enumerate(stat_names)
                           for group_name in group_names]
            print '\t'.join([''] + header_vals)
        col_vals = ['%.10f' % regressions[mi][gi, si + 2]
                    for si, stat_name in enumerate(stat_names)
                    for gi, group_name in enumerate(group_names)]
        print '\t'.join([measure_name] + col_vals)


def plot_stat_distributions(stats, group_names):
    # Show the distribution of stats, to see if there are
    # reliable group differences.
    fh4 = plt.figure(figsize=(18, 8))
    lbls = ['%s vs. %s' % (gn1, gn2)
            for gi, gn1 in enumerate(group_names)
            for gn2 in group_names[gi + 1:]]
    pi = 1
    for si, stat_name in enumerate(['mean', 'var']):
        stat_vals = np.asarray([ss[:, si * 2] for ss in stats])
        pvals = np.asarray([ss[:, si * 2 + 1] for ss in stats])
        for li, lbl in enumerate(lbls):
            ax1 = fh4.add_subplot(2, len(lbls), pi)
            ax1.hist(pvals[:, li], [0.0001, 0.001, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60])
            ax1.set_title(lbl)
            pi += 1
    equalize_xlims(fh4)
    equalize_ylims(fh4)


def loop_show_asymmetry(prefix,
                        grouping_prop_names=['Gender', 'FDH_23_Handedness_Prtcpnt'],
                        plots=['regressions', 'distributions']):
    """ Loop over all properties to show asymmetry."""
    data = load_PING_data()

    # Remove 'fuzzy':
    keys = [key for key in data.keys() if 'fuzzy' not in key]

    # Process & plot the data.
    stats = []
    regressions = []
    measure_names = get_twohemi_keys(keys, prefix=prefix)
    for pi, prop_name in enumerate(measure_names):
        # print "Comparing %d (%s)..." % (pi, prop_name)
        gn, ss, rv = compare_group_asymmetry(data, prop_name=prop_name, plots=plots,
                                         grouping_prop_names=grouping_prop_names)
        stats.append(ss)
        regressions.append(rv)

    if 'regression_stats' in plots:
        dump_regressions_csv(regressions,
                             group_names=gn,
                             measure_names=measure_names)

    if 'stat_distributions' in plots:
        plot_stat_distributions(stats, group_names=gn)

    plt.show()

if __name__ == '__main__':
    import warnings
    warnings.warn('Code to group by handedness or gender should be extracted and generalized.')
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

    loop_show_asymmetry(prefix=['MRI_cort_thick', 'MRI_cort_area', 'DTI_fiber_vol'],  # area_ctx_rh_frontalpole'],
                        grouping_prop_names=['FDH_23_Handedness_Prtcpnt'],
                        plots=['regression_stats'])
