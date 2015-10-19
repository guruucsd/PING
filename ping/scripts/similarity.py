"""
Similarity matrix comparisons PING data.
"""
import os
import simplejson
from collections import OrderedDict
from functools import partial

import numpy as np
from matplotlib import pyplot as plt

from ..ping.analysis.similarity import (compare_similarity_vectors,
                                        compute_similarity_vectors,
                                        visualize_similarity_matrices)
from ..ping.apps import PINGSession
from ..research.apps import ResearchArgParser
from ..research.asymmetry import is_ai_key
from ..research.data import get_all_data, dump_to_json
from ..research.plotting import show_plots


def do_similarity(prefixes, metric='partial-correlation', measures=None,
                  atlas='desikan', username=None, passwd=None,
                  output_format='matplotlib', data_dir='data', output_dir='data'):

    # Get prefix
    prefix_filter_fn = lambda k, v: np.any([k.startswith(p) for p in prefixes])

    # Load and filter the data
    p_data = get_all_data(atlas, username=username, passwd=passwd, data_dir=data_dir)
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
        measures = list(filt_fns.keys())

    # Get the key locations
    filt_fns = OrderedDict(((k, filt_fns[k]) for k in measures))

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
        class_label, label = [(p, p_data.get_nonhemi_key(key)[(len(p)):])
                              for p in prefixes
                              if key.startswith(p)][0]
        class_labels.append(class_label[:25])
        labels.append(p_data.get_anatomical_name(label))

    # Compare matrices (printed)
    compare_similarity_vectors(sim_dict)

    if output_format in ['json']:
        # Dump each prefix separately
        for ki, sim_mat in enumerate(sim_dict.values()):
            json_file = '%s.json' % (','.join(prefixes))
            json_file = os.path.join(output_dir, json_file)

            if sim_mat.ndim == 1:
                sim_mat = scipy.spatial.distance.squareform(sim_mat)
            sim_mat[np.eye(sim_mat.shape[0], dtype=bool)] = 0

            out_dict = dict()
            for li, lbl in enumerate(labels):
                out_dict[lbl] = dict(zip(labels, sim_mat[li]))

            dump_to_json(out_dict, json_file, klass=p_data.__class__)

    else:
        #  Display the similarity matrices.
        ax = visualize_similarity_matrices(sim_dict, labels=labels,
                                           class_labels=class_labels,
                                           output_format=output_format)
        ax.name = '%s' % ','.join(prefixes)
        show_plots(output_format, ax=ax, output_dir=output_dir)


if __name__ == '__main__':
    parser = ResearchArgParser(description="Compare asymmetry correlation "
                               "matrix with LH/RH structural covariance matrices.",
                               common_args=['prefixes',
                                            'atlas', 'username', 'passwd',
                                            'data-dir', 'output-dir'])
    parser.add_argument('metric', choices=['correlation', 'partial-correlation'],
                        nargs='?', default='partial-correlation')
    parser.add_argument('measures', choices=['Left Hemisphere',
                                            'Right Hemisphere',
                                            'Asymmetry Index',
                                            'all'],
                        nargs='?', default='all')
    parser.add_argument('--output-format', nargs='?', default='matplotlib',
                        choices=['matplotlib', 'mpld3', 'json',
                                 'bokeh', 'bokeh-silent'])
    args = parser.parse_args()
    args.measures = args.measures.split(',')
    do_similarity(**vars(args))
