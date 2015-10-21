"""
Similarity matrix comparisons PING data.
"""
import os
from collections import OrderedDict

import numpy as np

from ..research.apps import ResearchArgParser
from ..research.data import get_all_data, dump_to_json


def do_multivariate(prefixes, atlas='desikan', username=None, passwd=None,
                    data_dir='data', output_dir='.', output_format='json',
                    verbose=0, pc_thresh=0.05):
    # Select PCs
    filter_fns = dict([(p, lambda k, v: k.startswith(p))
                       for p in prefixes])
    data = get_all_data(atlas, filter_fns=filter_fns,
                        username=username, passwd=passwd,
                        data_dir=data_dir, verbose=verbose,
                        pc_thresh=pc_thresh)

    # Get the PC keys
    pc_keys = [k for k in data.data_dict.keys() if '_PC' in k]
    pcs = data.pcs

    # We can show PC loadings over subjects (mean, std),
    # as well as the components themselves.
    #
    # Brain 1: PC
    # Brain 2: mean, std (across subjects)
    #
    # Push this all into a single json.

    # First, compute summary stats from subject data and expand back onto each PC
    for pi, key in enumerate(pc_keys):
        pc_lbl = 'PC%d' % pi
        prefix = key[:-(len(pc_lbl)+1)]

        out_dict = dict()
        for alg in ['mean', 'std', 'pc']:
            subj_data = data.data_dict[key]
            components = np.asarray(pcs[pc_lbl].values())

            # Need to remove nan first
            projected_data = np.outer(subj_data, components)
            if alg == 'mean':
                summary_stat = np.nanmean(projected_data, axis=0)
            elif alg == 'std':
                summary_stat = np.nanstd(projected_data, axis=0)
            elif alg == 'pc':
                summary_stat = components

            # Stat map per ROI
            out_dict[alg] = dict(zip(pcs[pc_lbl].keys(), summary_stat))

        # now reverse mapping, so ROIs are outside, and mean/std/pc is inside.
        # { roi1: { mean: 1, std: 2, pc: 1}, roi2: ...}
        alg_keys = out_dict.keys()
        od_keys = out_dict.values()[0].keys()
        roi_keys = [k[len(prefix):-3] for k in od_keys]

        roi_vals = [dict(zip(alg_keys, [vals[k] for vals in out_dict.values()]))
                    for k in od_keys]

        json_file = os.path.join(output_dir, '%s.json' % (key))
        dump_to_json(dict(zip(roi_keys, roi_vals)), json_file, data.__class__)


if __name__ == '__main__':
    parser = ResearchArgParser(description="Compute PCA over asymmetries",
                               common_args=['prefixes',
                                            'atlas', 'username', 'passwd',
                                            'data-dir', 'output-dir',
                                            'verbose'])
    parser.add_argument('--output-format', nargs='?', default='json',
                        choices=['json', 'text'])
    parser.add_argument('--pc-thresh', type=float, default=0.05)
    args = parser.parse_args()
    do_multivariate(**vars(args))
