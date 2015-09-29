"""
Similarity matrix comparisons PING data.
"""
from argparse import ArgumentParser
from collections import OrderedDict
from functools import partial

import os
import numpy as np
from matplotlib import pyplot as plt

from ping.analysis.similarity import (compare_similarity_vectors,
                                      compute_similarity_vectors,
                                      visualize_similarity_matrices)
from ping.apps import PINGSession
from research.asymmetry import is_ai_key
from research.data import get_all_data


def do_similarity(prefix, metric='partial-correlation', measures=None,
                  dataset='ping', username=None, passwd=None,
                  plotengine='matplotlib'):

    # Get prefix
    prefix = prefix.split(',')
    prefix_filter_fn = lambda k, v: np.any([k.startswith(p) for p in prefix])

    # Load and filter the data
    p_data = get_all_data(dataset, username=username, passwd=passwd)
    p_data.filter(prefix_filter_fn)

    # Determine filters for selecting the similarity groupings.
    filt_fns = OrderedDict((
        ('Asymmetry Index', lambda key: (is_ai_key(key) and np.all([substr not in key for substr in ['_TOTAL', 'LH_PLUS_RH']])) or p_data.is_nonimaging_key(key)),
        ('Left Hemisphere', lambda key: p_data.which_hemi(key) == 'lh' or p_data.is_nonimaging_key(key)),
        ('Right Hemisphere', lambda key: p_data.which_hemi(key) == 'rh' or p_data.is_nonimaging_key(key)),
        ('all', lambda key: True)))
    filt_fns['LHAI'] = partial(lambda key, f1, f2: f1(key) or f2(key),
                               f1=filt_fns['Asymmetry Index'],
                               f2=filt_fns['Left Hemisphere'])

    # Get measures
    if measures is None:
        key_locations = list(filt_fns.keys())
    else:
        key_locations = measures.split(',')

    # Get the key locations
    filt_fns = OrderedDict(((k, filt_fns[k]) for k in key_locations))

    # Do the similarity computation; order by anatomy.
    sim_dict, sim_keys = compute_similarity_vectors(p_data.data_dict,
                                                    filt_fns=filt_fns,
                                                    sort_fn=p_data.anatomical_sort,
                                                    metric=metric)

    # Split the keys into the class (the prefix)
    # an anatomical label from the rest of the key.
    good_keys = next(iter(sim_keys.values()))
    labels = []
    class_labels = []
    for ki, key in enumerate(good_keys):
        class_label, label = [(p, p_data.get_nonhemi_key(key)[(len(p) + 1):])
                              for p in prefix
                              if key.startswith(p)][0]
        class_labels.append(class_label[:25])
        labels.append(p_data.get_anatomical_name(label)[:25])

    # Compare matrices (printed)
    compare_similarity_vectors(sim_dict)

    #  Display the similarity matrices.
    ax = visualize_similarity_matrices(sim_dict, labels=labels,
                                       class_labels=class_labels,
                                       plotengine=plotengine)
    # ax.get_figure().suptitle(', '.join([p_data.prefix2text(p) for p in prefix]), fontsize=24)

    # for key, mat in sim_dict.items():
    #     import scipy.spatial
    #     if mat.ndim == 1:
    #         mat = scipy.spatial.distance.squareform(mat)
    #     evals, evecs = np.linalg.eig(mat)
    #     # from research.multivariate import report_loadings
    #     # report_loadings(evals=evals, evecs=evecs, labels=np.asarray(labels))
    #     # Now print loadings, according to multivariate...
    if plotengine == 'mpld3':
        import mpld3
        mpld3.show()
    elif plotengine == 'matplotlib':
        plt.show()
    elif plotengine == 'bokeh':
        import bokeh.plotting
        bokeh.plotting.show(ax)


if __name__ == '__main__':
    parser = ArgumentParser(description="Compare asymmetry correlation"
                                        " matrix with LH/RH structural"
                                        " covariance matrices.\n")
    parser.add_argument('prefix', help="comma-separated list of prefixes to"
                                       " include in the analysis")
    parser.add_argument('metric', choices=['correlation', 'partial-correlation'],
                        nargs='?', default='partial-correlation')
    parser.add_argument('measures', choices=['Left Hemisphere',
                                            'Right Hemisphere',
                                            'Asymmetry Index',
                                            'all'],
                        nargs='?', default='all')
    parser.add_argument('--dataset', choices=['ping', 'destrieux'],
                        nargs='?', default='ping')
    parser.add_argument('--plotengine', choices=['matplotlib', 'mpld3', 'bokeh'],
                        nargs='?', default='matplotlib')
    parser.add_argument('--username', nargs='?',
                        default=PINGSession.env_username())
    parser.add_argument('--password', nargs='?',
                        default=PINGSession.env_passwd(),
                        dest='passwd')
    args = parser.parse_args()
    do_similarity(**vars(args))
