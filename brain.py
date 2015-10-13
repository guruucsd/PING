"""
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
from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
from flask import send_from_directory
from matplotlib import pyplot as plt
from six import string_types

import roygbiv
import roygbiv.server
from ping.analysis.similarity import is_bad_key
from ping.apps import PINGSession
from research.data import get_all_data, keytype2label
from research.plotting import show_plots
from scatter import compute_key_data

def do_roygbiv(prefix, key,
               dataset='desikan', username=None, passwd=None,
               output_format='json', subjects_dir=os.environ.get('SUBJECTS_DIR'),
               surface_type='pial', hemi='lh', subject='fsaverage',
               sample_rate=1., force=False):

    # Load the data (should group, but ... later.),
    # then filter by prefix
    data = get_all_data(dataset, username=username, passwd=passwd)
    data = data.filter(lambda k, v: np.any([k.startswith(p)
                                            for p in prefix.split(',')]))
    data = data.filter(lambda k, v: 'fuzzy' not in k)

    # Create the parcels
    atlas = dataset
    fsavg_path = os.path.join(subjects_dir, subject)
    surface_file = os.path.join(fsavg_path, 'surf', '%s.%s' % (hemi, surface_type))
    label_file = roygbiv.atlas2aparc(atlas, hemi=hemi)
    label_file = os.path.join(fsavg_path, 'label', label_file)
    output_dir = os.path.join('data', subject, atlas, surface_type)
    json_file = prefix + '%s_files_to_load.json' % hemi

    if force or not os.path.exists(os.path.join(output_dir, json_file)):
        roygbiv.freesurfer_annot_to_vtks(surface_file=surface_file,
                                         label_file=label_file,
                                         output_stem='%s%s_' % (prefix, hemi),
                                         json_file=json_file,
                                         output_dir=output_dir,
                                         sample_rate=sample_rate,
                                         force=force)

    def strip_prefix(key, prefix):
        if key.startswith(prefix):
            return key[len(prefix):]
        return key

    def map_colors(values):
        maxval = np.max(np.abs(values))
        minval = 0. if np.all(values >= 0) else -maxval
        norm = mpl.colors.Normalize(vmin=minval, vmax=maxval)

        cmap = cm.coolwarm
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        colors = [m.to_rgba(val)[0:3] for val in values]

        return colors

    if output_format in ['json', 'flask']:
        scatter_data1 = compute_key_data(data.data_dict, key=key)
        scatter_data2 = compute_key_data(data.data_dict, key='AI:std')
        scatter_data3 = compute_key_data(data.data_dict, key='LH_PLUS_RH:mean')
        non_prefixed_keys = [strip_prefix(key, prefix) for key in scatter_data1]
        anat_keys = [data.get_nonhemi_key(key) for key in non_prefixed_keys]
        labels = [data.get_anatomical_name(key) for key in anat_keys]
        values = zip(scatter_data1.values(), scatter_data2.values(), scatter_data3.values())
        colors = map_colors(np.asarray(scatter_data1.values()))
        out_dict = dict(names=dict(zip(anat_keys, labels)),
                        values=dict(zip(anat_keys, values)),
                        colors=dict(zip(anat_keys, colors)))
        print json_file
        roygbiv.add_metadata(out_dict, json_file=json_file, output_dir=output_dir)

        if output_format == 'flask':
            data_dir = os.path.abspath('data')
            rgb_dir = os.path.dirname(os.path.abspath(roygbiv.__file__))
            web_dir = os.path.join(rgb_dir, '..', 'web')
            app = roygbiv.server.make_server(data_dir=data_dir)

            @app.route('/<path:dataset>/<path:atlas>/<path:surface>/<path:prefix>/data/<path:path>')
            def send_data_specific_new(dataset, atlas, surface, prefix, path):
                if path.endswith('json'):
                    fn = os.path.basename(path)
                    path = os.path.join(os.path.dirname(path), prefix + fn)
                cur_dir = os.path.join(data_dir, dataset, atlas, surface)
                print cur_dir, path
                return send_from_directory(cur_dir, path)

            @app.route('/<path:dataset>/<path:atlas>/<path:surface>/<path:prefix>/<path:html_file>')
            def send_allspecific_new(dataset, atlas, surface, prefix, html_file):
                if html_file == '':
                    html_file = 'index'
                print web_dir, html_file
                return send_from_directory(web_dir, html_file)

            app.run()

    else:
        raise NotImplementedError()



if __name__ == '__main__':
    axis_choices = ['AI:mean', 'AI:std',
                    'LH_PLUS_RH:mean', 'LH_PLUS_RH:std']
    parser = ArgumentParser(description="Scatter plot on any two data"
                            " arrays, with additional data arrays that"
                            " optionally control marker size and color.")
    parser.add_argument('prefix', help="comma-separated list of prefixes to"
                                       " include in the analysis")
    parser.add_argument('key', choices=axis_choices)
    parser.add_argument('--dataset', choices=['desikan', 'destrieux'],
                        nargs='?', default='desikan')
    parser.add_argument('--output-format', choices=['flask', 'json'],
                        nargs='?', default='json')
    parser.add_argument('--hemi', choices=['lh', 'rh'],
                        nargs='?', default='lh')
    parser.add_argument('--sample-rate', nargs='?', default=1.)
    parser.add_argument('--subject', nargs='?', default='fsaverage')
    parser.add_argument('--surface-type', choices=['pial', 'inflated'],
                        nargs='?', default='pial')
    parser.add_argument('--username', nargs='?',
                        default=PINGSession.env_username())
    parser.add_argument('--password', nargs='?',
                        default=PINGSession.env_passwd(),
                        dest='passwd')
    parser.add_argument('--force', nargs='?',
                        default=False, choices=[False, True])
    args = parser.parse_args()
    do_roygbiv(**vars(args))
