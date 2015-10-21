"""
File for investigating asymmetry from PING data, based on each subject's
asymmetry index
"""
import sys
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

from ..ping.apps import PINGSession
from ..ping.data import PINGData
from ..ping.utils import do_and_plot_regression
from ..ping.utils.plotting import (plot_symmetric_matrix_as_triangle,
                                   equalize_xlims, equalize_ylims,
                                   plot_normalized_hist)
from ..research.apps import ResearchArgParser
from ..research.asymmetry import get_asymmetry_index
from ..research.data import get_all_data
from ..research.plotting import show_plots


def get_groupings(data, grouping_keys):
    """ Group data by the unique values of `grouping_key`"""

    if not isinstance(grouping_keys, list):
        grouping_keys = [grouping_keys]

    # Group
    grouping_index = group_names = []
    for gpn in grouping_keys:
        grouping_data = np.asarray(data[gpn].tolist())
        if grouping_data.dtype == float:
            grouping_data = np.round(grouping_data).astype(int)
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


def compute_group_asymmetry(data, xaxis_key, yaxis_key, grouping_keys):
    group_names, grouping_index = get_groupings(data.data_dict, grouping_keys)
    x_data = np.asarray(data.data_dict[xaxis_key].tolist())
    prop_ai = data.data_dict[yaxis_key]
    # prop_ai = get_asymmetry_index(data, yaxis_key)

    group_x = []
    group_y = []
    for gi, group_name in enumerate(group_names):

        # Index the current group
        if group_name == 'all':
            idx = np.ones((len(data.data_dict.values()[0]),), dtype=bool)
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

        group_x.append(cur_x)
        group_y.append(group_ai)
    return group_names, group_x, group_y


def plot_regressions(group_names, group_x, group_y, plotengine='matplotlib',
                     xaxis_key=None, yaxis_key=None):

    n_subplots = len(group_names)
    n_rows = 1  # int(np.round(np.sqrt(n_subplots)))
    n_cols = n_subplots  # int(np.ceil(n_subplots / float(n_rows)))

    fh1 = plt.figure(figsize=(18, 7))
    # fh1.suptitle('%s regression' % yaxis_key)

    regressions = []
    for gi, (group_name, gx, gy) in enumerate(zip(group_names, group_x, group_y)):
        regressions.append(scipy.stats.linregress(gx, gy))

        # Plot the regression result
        params = dict(xlabel=xaxis_key, ylabel='Asymmetry Index (LH - RH)',
                      title='Group: %s (n=%d)' % (group_name, len(gy)))
        if gi > 0:
            del params['ylabel']
        ax1 = fh1.add_subplot(n_rows, n_cols, gi + 1)
        ax1 = fh1.gca()
        do_and_plot_regression(gx, gy, ax=ax1, colori=gi,
                               show_std=(len(gx) > 200),
                               plotengine=plotengine, **params)
        #ax1.set_title(measure_key)  # ax1.get_title().split('\n')[0])
    regressions = np.asarray(regressions)
    ax1.legend(group_names)
    return regressions


def plot_distributions(group_names, group_x, group_y, xaxis_key, yaxis_key, n_bins=15):
    n_subplots = len(group_names)
    n_rows = 1  # int(np.round(np.sqrt(n_subplots)))
    n_cols = n_subplots  # int(np.ceil(n_subplots / float(n_rows)))

    fh2 = plt.figure(figsize=(18, 7))
    fh2.suptitle('%s distributions' % yaxis_key)
    bins = np.linspace(np.asarray([ys.min() for ys in group_y]).min(),
                       np.asarray([ys.max() for ys in group_y]).max(),
                       n_bins + 2)[1:-1]

    for gi, (group_name, gx, gy) in enumerate(zip(group_names, group_x, group_y)):
        # Plot the regression result
        params = dict(xlabel=xaxis_key, ylabel='Asymmetry Index (LH - RH)',
                      title='Group: %s (n=%d)' % (group_name, len(gy)))

        # Plot the regression result
        ax2 = fh2.add_subplot(n_rows, n_cols, gi + 1)
        plot_normalized_hist(gy, ax2, bins=bins)
        ax2.set_title(params['title'])
        ax2.set_xlabel(params['ylabel'])
        ax2.set_ylim([0, 0.25])
    equalize_xlims(fh2)
    equalize_ylims(fh2)


def do_stats(group_names, group_x, group_y):
    # stats[:, n]: 0:mean_stat: 1:mean_pval; 2:std_stat, 3:std_pval
    # shape: n_compares x 4
    stats = np.asarray([(scipy.stats.ttest_ind(gsamps1, gsamps2) +
                         scipy.stats.levene(gsamps1, gsamps2))
                        for gi, gsamps1 in enumerate(group_y)
                        for gsamps2 in group_y[gi + 1:]])

    # Test whether variances differ
    dist_mat = scipy.spatial.distance.squareform(stats[:, 2])
    sig_mat = stats[:, 3] <= (0.05 / stats.shape[0])

    fh3 = plt.figure()
    fh3.suptitle(str(['%.2e' % s for s in stats[:, 3]]))

    ax1 = fh3.add_subplot(1, 2, 1)
    plot_symmetric_matrix_as_triangle(dist_mat, ax=ax1, labels=group_names)

    ax2 = fh3.add_subplot(1, 2, 2, axisbg=fh3.get_facecolor())
    plot_symmetric_matrix_as_triangle(sig_mat, ax=ax2, labels=group_names)



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


def do_grouping(prefixes, grouping_keys, xaxis_key='Age_At_IMGExam',
                plots='regressions',
                atlas='desikan', username=None, passwd=None,
                data_dir='data', output_dir='.', output_type='matplotlib'):
    """ Loop over all properties to show asymmetry."""

    data = get_all_data(atlas, username=username, passwd=passwd, data_dir=data_dir)
    data.filter(lambda k, v: 'fuzzy' not in k)  # Remove 'fuzzy'
    data.filter([partial(lambda k, v, p: (k.startswith(p) or
                                          k in grouping_keys or
                                          k == xaxis_key),
                         p=p)
                 for p in prefixes])

    # Process & plot the data.
    measure_keys = [(k) for k in data.data_dict.keys()
                        if k.startswith(p) and k.endswith('_AI')]
    for pi, yaxis_key in enumerate(sorted(measure_keys)):
        print("Comparing %d (%s)..." % (pi, yaxis_key))

        kwargs = dict(xaxis_key=xaxis_key, yaxis_key=yaxis_key)
        group_names, group_x, group_y = compute_group_asymmetry(data, grouping_keys=grouping_keys, **kwargs)

        kwargs.update(dict(group_names=group_names, group_x=group_x, group_y=group_y))
        if 'regressions' in plots:
            plot_regressions(**kwargs)
        if 'distributions' in plots:
            plot_distributions(**kwargs)
        if 'stats' in plots:
            do_stats(group_names, group_x, group_y)

    # Outside of loop
    if 'regression_stats' in plots:
        dump_regressions_csv(regressions,
                             group_names=group_names,
                             measure_names=measure_keys)

        plot_regressions_scatter(regressions,
                                 group_names=group_names,
                                 measure_names=measure_keys)

    if 'stat_distributions' in plots:
        plot_stat_distributions(stats, group_names=group_names)

    show_plots(plotengine=output_type, output_dir=output_dir)


if __name__ == '__main__':
    parser = ResearchArgParser(description="Produce plots for each group.",
                               common_args=['prefixes',
                                            'atlas', 'username', 'passwd',
                                            'data-dir', 'output-dir'])
    parser.add_argument('grouping_keys', choices=['Gender', 'FDH_23_Handedness_Prtcpnt'])
    parser.add_argument('--xaxis_key', help="spreadsheet value to regress against.",
                        nargs='?', default='Age_At_IMGExam')
    parser.add_argument('--plots', choices=['regressions', 'distributions',
                                          'stats', 'regression_stats',
                                          'stat_distributions'],
                        nargs='?', default='regressions',
                        help="comma-separated list of plots")
    parser.add_argument('--output-type', choices=['matplotlib', 'mpld3'],
                        nargs='?', default='matplotlib')

    args = parser.parse_args()
    args.grouping_keys = args.grouping_keys.split(',')
    args.plots = args.plots.split(',')

    do_grouping(**vars(args))
