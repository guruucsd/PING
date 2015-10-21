i"""
Various scatter plots

Goal is to have:
2D: typical scatter
3D: z is the size of the marker
4D: add in the color of the marker.

Should take ordered parameters for data keys on the input,
function should take keyword args.
"""

import os
import simplejson

import numpy as np
from flask import send_from_directory
from matplotlib import pyplot as plt
from six import string_types

import roygbiv
import roygbiv.server
from .scatter import compute_key_data
from ..ping.analysis.similarity import is_bad_key
from ..ping.apps import PINGSession
from ..research.apps import ResearchArgParser
from ..research.data import get_all_data, keytype2label, map_colors, strip_prefix
from ..research.plotting import show_plots


def do_roygbiv(prefix, key,
               atlas='desikan', username=None, passwd=None,
               output_format='json', subjects_dir=os.environ.get('SUBJECTS_DIR'),
               surface_type='pial', hemi='lh', subject='fsaverage',
               sample_rate=1., force=False, data_dir='data', output_dir='data'):

    # Load the data (should group, but ... later.),
    # then filter by prefix
    data = get_all_data(atlas, username=username, passwd=passwd, data_dir=data_dir)
    data = data.filter(lambda k, v: np.any([k.startswith(p)
                                            for p in prefix.split(',')]))
    data = data.filter(lambda k, v: 'fuzzy' not in k)

    # Create the parcels
    fsavg_path = os.path.join(subjects_dir, subject)
    surface_file = os.path.join(fsavg_path, 'surf', '%s.%s' % (hemi, surface_type))
    label_file = roygbiv.atlas2aparc(atlas, hemi=hemi)
    label_file = os.path.join(fsavg_path, 'label', label_file)
    output_dir = os.path.join(output_dir, subject, atlas, surface_type)
    json_file = prefix + '%s_files_to_load.json' % hemi

    if force or not os.path.exists(os.path.join(output_dir, json_file)):
        roygbiv.freesurfer_annot_to_vtks(surface_file=surface_file,
                                         label_file=label_file,
                                         output_stem='%s_' % (hemi),
                                         json_file=json_file,
                                         output_dir=output_dir,
                                         sample_rate=sample_rate,
                                         force=force)

    if output_format in ['json', 'flask']:
        assert key == 'AI:mean'
        scatter_data1 = compute_key_data(data.data_dict, key='AI:mean')
        scatter_data2 = compute_key_data(data.data_dict, key='AI:std')
        scatter_data3 = compute_key_data(data.data_dict, key='LH_PLUS_RH:mean')

        if hemi == 'rh':
            for key in scatter_data1:
                scatter_data1[key] = -scatter_data1[key]
                scatter_data2[key] = -scatter_data2[key]


        # Compute keys / labels
        non_prefixed_keys = [strip_prefix(key, prefix) for key in scatter_data1]
        anat_keys = [data.get_nonhemi_key(key) for key in non_prefixed_keys]
        labels = [data.get_anatomical_name(key) for key in anat_keys]

        # Compute values
        values = zip(scatter_data1.values(), scatter_data2.values(), scatter_data3.values())
        colors = map_colors(np.asarray(scatter_data1.values()))

        # Aggregate
        out_dict = dict(names=dict(zip(anat_keys, labels)),
                        values=dict(zip(anat_keys, values)),
                        colors=dict(zip(anat_keys, colors)))
        print json_file
        roygbiv.add_metadata(out_dict, json_file=json_file, output_dir=output_dir)

        if output_format == 'flask':
            web_dir = os.path.dirname(os.path.abspath(roygbiv.web))
            app = roygbiv.server.make_server(data_dir=output_dir)

            @app.route('/data/<path:path>')
            def send_data_specific_new(path):
                return send_from_directory(cur_dir, path)

            @app.route('/<path:dataset>/<path:atlas>/<path:surface>/<path:html_file>')
            def send_allspecific_new(dataset, atlas, surface, prefix, html_file):
                if html_file == '':
                    html_file = 'index'
                return send_from_directory(web_dir, html_file)

            app.run()

    else:
        raise NotImplementedError()


if __name__ == '__main__':
    parser = ResearchArgParser(description="Scatter plot on any two data"
                               " arrays, with additional data arrays that"
                               " optionally control marker size and color.",
                               common_args=['prefix', 'key',
                                            'atlas', 'hemi',
                                            'force', 'data-dir', 'output-dir'])
    parser.add_argument('--sample-rate', nargs='?', default=1.)
    parser.add_argument('--subject', nargs='?', default='fsaverage')
    parser.add_argument('--surface-type', choices=['pial', 'inflated'],
                        nargs='?', default='pial')
    parser.add_argument('--output-format', choices=['json', 'flask'],
                        nargs='?', default='json')
    args = parser.parse_args()
    do_roygbiv(**vars(args))
