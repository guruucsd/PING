"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import os
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from ping.data import PINGData
from ping.utils import do_and_plot_regression
from ping.utils.plotting import (plot_symmetric_matrix_as_triangle,
                                 equalize_xlims, equalize_ylims,
                                 plot_normalized_hist)
from research.asymmetry import get_asymmetry_index


def make_groups(data, grouping_keys):
    """ Group data by the unique values of `grouping_key`"""

    if not isinstance(grouping_keys, list):
        grouping_keys = [grouping_keys]

    # Group
    grouping_index = group_names = []
    for gpn in grouping_keys:
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


def compare_group_asymmetry(data, key, grouping_keys, plots):
    """ Groups data according to grouping_keys, computes
    asymmetry index for key, and graphs."""

    group_names, grouping_index = make_groups(data, grouping_keys)

    age_data = np.asarray(data['Age_At_IMGExam'].tolist())
    prop_ai = get_asymmetry_index(data, key)

    n_subplots = len(group_names)
    n_rows = int(np.round(np.sqrt(n_subplots)))
    n_cols = int(np.ceil(n_subplots / float(n_rows)))

    if 'regressions' in plots:
        fh1 = plt.figure(figsize=(18, 10))
        fh1.suptitle('%s regression' % key)
    if 'distributions' in plots:
        fh2 = plt.figure(figsize=(18, 10))
        fh2.suptitle('%s distributions' % key)
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
            plot_normalized_hist(group_ai, ax2, bins=bins)
            ax2.set_title(params['title'])
            ax2.set_xlabel(params['ylabel'])
            ax2.set_ylims([0, 0.25])
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
        sig_mat = stats[:, 3] <= (0.05 / stats.shape[0])

        fh3 = plt.figure()
        fh3.suptitle(str(['%.2e' % s for s in stats[:, 3]]))

        ax1 = fh3.add_subplot(1, 2, 1)
        plot_symmetric_matrix_as_triangle(dist_mat, ax=ax1, labels=group_names)

        ax2 = fh3.add_subplot(1, 2, 2, axisbg=fh3.get_facecolor())
        plot_symmetric_matrix_as_triangle(sig_mat, ax=ax2, labels=group_names)

    if 'distributions' in plots:
        equalize_xlims(fh2)
        equalize_ylims(fh2)

    return group_names, stats, regressions, group_samples


def dump_regressions_csv(regressions, group_names, measure_names):
    # Dump a tsv of group rvals, pvals, and coeff of variation
    n_measures = len(regressions)
    stat_names = ['rval', 'pval', 'coeff_of_var']
    for mi, measure_name in enumerate(sorted(measure_names)):
        if mi == 0:
            header_vals = ['%s_%s' % (group_name, stat_name)
                           for si, stat_name in enumerate(stat_names)
                           for group_name in group_names]
            print('\t'.join([''] + header_vals))

        row_vals = ['%.10f' % regressions[mi][gi, si + 2]
                    for si, stat_name in enumerate(stat_names)
                    for gi, group_name in enumerate(group_names)]
        print('\t'.join([measure_name] + row_vals))


def plot_regressions_scatter(regressions, group_names, measure_names):
    # Dump a tsv of group rvals, pvals, and coeff of variation
    n_measures = len(regressions)
    stat_names = ['rval', 'pval']
    color_arr = ['b', 'g', 'r', 'k', 'y', 'c'][:len(group_names)]

    prefix = np.unique([m[:12] for m in measure_names])
    all_xvals = dict([(gn, []) for gn in group_names])
    all_yvals = dict([(gn, []) for gn in group_names])

    for p in prefix:
        prefix_idx = np.asarray([m.startswith(p) for m in measure_names])

        p_xvals = dict()
        p_yvals = dict()
        for gi, group_name in enumerate(group_names):
            # Plot rval vs. pval
            group_xvals = np.asarray([regressions[mi][gi, 2]
                                      for mi in np.nonzero(prefix_idx)[0]])
            group_yvals = np.asarray([regressions[mi][gi, 3]
                                      for mi in np.nonzero(prefix_idx)[0]])

            good_idx = ~np.logical_or(np.isnan(group_xvals),
                                      np.isnan(group_yvals))
            p_xvals[group_name] = group_xvals[good_idx]
            p_yvals[group_name] = group_yvals[good_idx]

            all_xvals[group_name] += group_xvals[good_idx].tolist()
            all_yvals[group_name] += group_yvals[good_idx].tolist()

        fh = plt.figure(figsize=(18, 8))
        ax1 = fh.add_subplot(1, 2, 1)
        plot_normalized_hist(
            [p_xvals[gn] for gn in group_names]).T,
            ax=ax1,
            bins=20,
            color=color_arr)
        ax1.legend(group_names)
        ax2 = fh.add_subplot(1, 2, 2)
        plot_normalized_hist(
            np.asarray([p_yvals[gn] for gn in group_names]).T,
            ax=ax2,
            bins=20,
            color=color_arr)
        fh.suptitle(p)

    # Bar plot (total)
    fh = plt.figure(figsize=(18, 8))
    ax1 = fh.add_subplot(1, 2, 1)
    plot_normalized_hist(np.asarray([all_xvals[gn] for gn in group_names]).T,
                         ax=ax1,
                         bins=50,
                         color=color_arr)
    ax1.legend(group_names)
    ax2 = fh.add_subplot(1, 2, 2)
    plot_normalized_hist(np.asarray([all_yvals[gn] for gn in group_names]).T,
                         ax=ax2,
                         bins=50,
                         color=color_arr)
    fh.suptitle('Over all measures')

    #ax = plt.figure(figsize=(18, 8)).gca()
    ax1.scatter(np.asarray([all_xvals[gn] for gn in group_names]).T,
                30 * np.asarray([all_yvals[gn] for gn in group_names]).T,
                c=color_arr)


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
            plot_normalized_hist(pvals[:, li], ax=ax1, bins=[0.0001, 0.001, 0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60])
            ax1.set_title(lbl)
            if li == 0:
                ax1.set_ylabel(stat_name)
            pi += 1
    equalize_xlims(fh4)
    equalize_ylims(fh4)

    means = np.asarray([ss[:, 0] for ss in stats])
    stds = np.asarray([ss[:, 2] for ss in stats])
    good_idx = ~np.logical_or(np.isnan(means.sum(1)), np.isnan(stds.sum(1)))

    import pdb; pdb.set_trace()
    fh = plt.figure(figsize=(18, 8))
    ax1 = fh.add_subplot(1, 2, 1)
    ax1.hist(means[good_idx], bins=25)
    ax1.set_xlabel('mean')

    ax2 = fh.add_subplot(1, 2, 2)
    ax2.hist(stds[good_idx], bins=25)
    ax2.set_xlabel('standard deviation')
    ax2.legend(group_names)

def loop_show_asymmetry(prefix,
                        grouping_keys=['Gender', 'FDH_23_Handedness_Prtcpnt'],
                        plots=['regressions', 'distributions']):
    """ Loop over all properties to show asymmetry."""
    data = PINGData()
    data.filter(lambda k, v: 'fuzzy' not in k)  # Remove 'fuzzy'
    data.filter([partial(lambda k, v, p: (k.startswith(p) or
                                          k in grouping_keys or
                                          k == 'Age_At_IMGExam'),
                         p=p)
                 for p in prefix])

    # Process & plot the data.
    stats = []
    regressions = []
    measure_keys = data.get_twohemi_keys()
    for pi, key in enumerate(sorted(measure_keys)):
        # print("Comparing %d (%s)..." % (pi, key))
        print(key)
        gn, ss, rv, gs = compare_group_asymmetry(data.data_dict, key=key, plots=plots,
                                                 grouping_keys=grouping_keys)
        stats.append(ss)
        regressions.append(rv)
        group_samples.append(gs)

    if 'regression_stats' in plots:
        dump_regressions_csv(regressions,
                             group_names=gn,
                             measure_names=measure_keys)

        plot_regressions_scatter(regressions,
                                 group_names=gn, 
                                 measure_names=measure_keys)

    if 'stat_distributions' in plots:
        plot_stat_distributions(stats, group_names=gn)

    plt.show()

if __name__ == '__main__':
    import warnings
    warnings.warn('Code to group by handedness or gender should be extracted and generalized.')

    loop_show_asymmetry(prefix=PINGData.IMAGING_PREFIX, #['MRI_cort_thick'],  # area_ctx_rh_frontalpole'],
                        grouping_keys=['Gender'],
                        plots=['stat_distributions'])
