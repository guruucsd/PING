"""
Similarity matrix comparisons PING data.
"""
from collections import OrderedDict
from functools import partial

import numpy as np
from matplotlib import pyplot as plt

from ping.analysis.similarity import (compare_similarity_vectors,
                                      compute_similarity_vectors,
                                      visualize_similarity_matrices)
from research.asymmetry import is_ai_key
from research.data import get_all_data


def do_usage(args, error_msg=None):
    if error_msg is not None:
        print("*** ERROR *** : %s" % error_msg)
    print("\nUsage: %s prefix metric which_matrices" % __file__)
    print("\tCompare asymmetry correlation matrix with LH/RH structural covariance matrices.")
    print("\n\tprefix: comma-separated list of prefixes to include in the analysis.")
    print("\tmetric: correlation or partial-correlation.")
    print("\twhat: Left Hemisphere, Right Hemisphere, Asymmetry Index, All")


def do_similarity(*args):
    if len(args) > 3:
        do_usage(args, error_msg="Too many arguments.")
        return

    elif len(args) < 1:
        do_usage(args, error_msg="Too few arguments.")
        return

    # Get prefix
    prefix = args[0].split(',')
    prefix_filter_fn = lambda k, v: np.any([k.startswith(p) for p in prefix])

    # Load and filter the data
    p_data = get_all_data()
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

    # Get metric
    metric = 'partial-correlation' if len(args) <= 1 else args[1]
    key_locations = list(filt_fns.keys()) if len(args) <= 2 else args[2].split(',')

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
        class_labels.append(class_label)
        labels.append(p_data.get_anatomical_name(label))

    # Compare matrices (printed)
    compare_similarity_vectors(sim_dict)

    #  Display the similarity matrices.
    ax = visualize_similarity_matrices(sim_dict, labels=labels,
                                       class_labels=class_labels)
    # ax.get_figure().suptitle(', '.join([p_data.prefix2text(p) for p in prefix]), fontsize=24)

    for key, mat in sim_dict.items():
        import scipy.spatial
        if mat.ndim == 1:
            mat = scipy.spatial.distance.squareform(mat)
        evals, evecs = np.linalg.eig(mat)
        # from research.multivariate import report_loadings
        # report_loadings(evals=evals, evecs=evecs, labels=np.asarray(labels))
        # Now print loadings, according to multivariate...
    plt.show()


if __name__ == '__main__':
    import sys
    do_similarity(*sys.argv[1:])
