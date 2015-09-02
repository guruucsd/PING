"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import os
import sys
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from ping.utils import do_and_plot_regression
from ping.utils.plotting import (plot_symmetric_matrix_as_triangle,
                                 equalize_xlims, equalize_ylims,
                                 plot_normalized_hist)
from research.asymmetry import get_asymmetry_index
from research.data import get_all_data
from research.grouping import get_groupings


def compare_group_asymmetry(data, xaxis_key, yaxis_key, grouping_keys, plots, measure_key):
    """ Groups data according to grouping_keys, computes
    asymmetry index for key, and graphs."""

    group_names, grouping_index = get_groupings(data, grouping_keys)

    x_data = np.asarray(data[xaxis_key].tolist())
    prop_ai = get_asymmetry_index(data, yaxis_key)

    n_subplots = len(group_names)
    n_rows = 1  # int(np.round(np.sqrt(n_subplots)))
    n_cols = n_subplots  # int(np.ceil(n_subplots / float(n_rows)))

    if 'regressions' in plots:
        fh1 = plt.figure(figsize=(18, 7))
        # fh1.suptitle('%s regression' % yaxis_key)
    if 'distributions' in plots:
        fh2 = plt.figure(figsize=(18, 7))
        fh2.suptitle('%s distributions' % yaxis_key)
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
        cur_x = x_data[idx]
        group_ai = prop_ai[idx]

        # Remove bad data
        bad_idx = np.logical_or(np.isnan(cur_x), np.isnan(group_ai))
        good_idx = np.logical_not(bad_idx)
        cur_x = cur_x[good_idx]
        group_ai = group_ai[good_idx]

        group_samples.append(group_ai)
        regressions.append(scipy.stats.linregress(cur_x, group_ai))

        # Plot the regression result
        params = dict(xlabel=xaxis_key, ylabel='Asymmetry Index (LH - RH)',
                      title='Group: %s (n=%d)' % (group_name, len(group_ai)))

        if 'regressions' in plots:
            if gi > 0:
                del params['ylabel']
            # ax1 = fh1.add_subplot(n_rows, n_cols, gi + 1)
            ax1 = fh1.gca()
            do_and_plot_regression(cur_x, group_ai, ax=ax1, colori=gi,
                                   show_std=(len(cur_x) > 200), **params)
            ax1.set_title(measure_key)  # ax1.get_title().split('\n')[0])

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
            np.asarray([p_xvals[gn] for gn in group_names]).T,
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

    fh = plt.figure(figsize=(18, 8))
    ax1 = fh.add_subplot(1, 2, 1)
    plot_normalized_hist(means[good_idx], ax=ax1, bins=25)
    ax1.set_xlabel('mean')

    ax2 = fh.add_subplot(1, 2, 2)
    plot_normalized_hist(stds[good_idx], ax=ax2, bins=25)
    ax2.set_xlabel('standard deviation')
    ax2.legend(group_names)


def loop_show_asymmetry(prefix,
                        grouping_keys=['Gender', 'FDH_23_Handedness_Prtcpnt'],
                        xaxis_key='Age_At_IMGExam',
                        plots='regressions'):
    """ Loop over all properties to show asymmetry."""
    data = get_all_data()
    data.filter(lambda k, v: 'fuzzy' not in k)  # Remove 'fuzzy'
    data.filter([partial(lambda k, v, p: (k.startswith(p) or
                                          k in grouping_keys or
                                          k == xaxis_key),
                         p=p)
                 for p in prefix])

    # Process & plot the data.
    stats = []
    regressions = []
    group_samples = []
    measure_keys = data.get_twohemi_keys()
    for pi, key in enumerate(sorted(measure_keys)):
        print("Comparing %d (%s)..." % (pi, key))
        gn, ss, rv, gs = compare_group_asymmetry(data.data_dict, xaxis_key=xaxis_key,
                                                 yaxis_key=key, plots=plots,
                                                 grouping_keys=grouping_keys,
                                                 measure_key=key)
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


def do_usage(args, error_msg=None):
    if error_msg is not None:
        print("*** ERROR *** : %s" % error_msg)
    print("\nUsage: %s prefixes group_keys [xaxis] [plots]" % __file__)
    print("\tProduce plots for each group")
    print("\n\tprefixes: simple selector for groups of measures to include. Popular choices include:")
    print("\t\tMRI_cort_area.ctx: ")
    print("\t\tMRI_cort_thick.ctx: ")
    print("\n\tgroups: comma separated list of keys for grouping. Popular ones include:")
    print("\t\tGender: ")
    print("\t\tFDH_23_Handedness_Prtcpnt: ")
    print("\txaxis: (optional) value to regress against (default=Age_At_IMGExam)")
    print("\tplots: (optional) list of plots (default=regressions); select from:")
    print("\t\tregressions:")
    print("\t\tdistributions:")
    print("\t\tstats:")
    print("\t\tregression_stats:")
    print("\t\tstat_distributions:")


def do_grouping(*args):
    if len(args) >= 5:
        do_usage(args, "Too many arguments.")
        return

    elif len(args) < 2:
        do_usage(args, "Too few arguments.")
        return

    prefix = args[0].split(',')
    grouping_keys = args[1].split(',')
    xaxis_key = args[2] if len(args) >= 3 else 'Age_At_IMGExam'
    plots = args[3].split(',') if len(args) >= 4 else 'regressions'

    loop_show_asymmetry(prefix=prefix,
                        grouping_keys=grouping_keys,
                        xaxis_key=xaxis_key,
                        plots=plots)


if __name__ == '__main__':
    import sys
    do_grouping(*sys.argv[1:])


